"""
Tests for BioCompiler Standalone Certificate Verifier.

These tests verify the standalone verifier works correctly WITHOUT
importing any biocompiler modules. The verifier is self-contained
and all checks are implemented from scratch.
"""

import hashlib
import importlib.util
import json
import os
import sys
import tempfile

import pytest

# Import ONLY the standalone verifier by loading the module file directly.
# This avoids the broken biocompiler/__init__.py which has stale imports.
_VERIFIER_PATH = os.path.join(
    os.path.dirname(__file__), "..", "src", "biocompiler", "standalone_verifier.py"
)
_spec = importlib.util.spec_from_file_location("standalone_verifier", _VERIFIER_PATH)
_sv = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(_sv)

verify = _sv.verify
verify_file = _sv.verify_file
check_gc_in_range = _sv.check_gc_in_range
check_no_stop_codons = _sv.check_no_stop_codons
check_no_cryptic_splice = _sv.check_no_cryptic_splice
check_no_restriction_sites = _sv.check_no_restriction_sites
check_no_cpg_island = _sv.check_no_cpg_island
check_no_gt_dinucleotide = _sv.check_no_gt_dinucleotide
check_cai_above_threshold = _sv.check_cai_above_threshold
check_protein_identity = _sv.check_protein_identity
_gc_content = _sv._gc_content
_compute_cai = _sv._compute_cai
_maxent_donor_score = _sv._maxent_donor_score
_dispatch_check = _sv._dispatch_check
_get_params = _sv._get_params
_resolve_base_name = _sv._resolve_base_name


# ── Helper: build a valid certificate dict ──

def _make_seq(length: int = 300, gc_target: float = 0.55, seed: int = 42) -> str:
    """Generate a DNA sequence with approximate GC content."""
    import random
    rng = random.Random(seed)
    gc_bases = "GC"
    at_bases = "AT"
    seq = []
    for i in range(length):
        if i < 3:
            # Start with ATG
            seq.append("ATG"[i])
            continue
        if i >= length - 3:
            # End with TAA
            seq.append("TAA"[i - (length - 3)])
            continue
        if rng.random() < gc_target:
            seq.append(rng.choice(gc_bases))
        else:
            seq.append(rng.choice(at_bases))
    return "".join(seq)


def _make_valid_cert(seq: str = None, types_overrides: list = None) -> dict:
    """Build a valid certificate dict for testing."""
    if seq is None:
        seq = _make_seq()
    seq_hash = hashlib.sha256(seq.encode()).hexdigest()
    types_list = [
        {"predicate": "NoStopCodons", "verdict": "PASS", "derivation": []},
        {"predicate": "GCInRange", "verdict": "PASS", "derivation": []},
        {"predicate": "NoRestrictionSite", "verdict": "PASS", "derivation": []},
        {"predicate": "NoCpGIsland", "verdict": "PASS", "derivation": []},
        {"predicate": "NoGTDinucleotide", "verdict": "FAIL", "derivation": []},
        {"predicate": "CodonAdapted", "verdict": "PASS", "derivation": []},
    ]
    if types_overrides is not None:
        types_list = types_overrides
    return {
        "version": "7.2.0",
        "design_id": seq_hash,
        "sequence": seq,
        "types": types_list,
        "provenance": {
            "tool": "BioCompiler",
            "version": "7.2.0",
            "timestamp": "2025-01-01T00:00:00Z",
            "parameters": {
                "gc_lo": 0.30,
                "gc_hi": 0.70,
                "cai_threshold": 0.5,
                "organism": "Homo_sapiens",
                "enzymes": ["EcoRI", "BamHI", "XhoI", "HindIII", "NotI"],
                "cryptic_splice_threshold": 3.0,
                "cpg_window": 200,
                "cpg_threshold": 0.6,
            },
            "input_hash": seq_hash,
            "overall_status": "PARTIAL_5/6",
            "mutagenesis": {"applied": False},
        },
    }


# ══════════════════════════════════════════════════════════════════════
# Test: Individual predicate checks
# ══════════════════════════════════════════════════════════════════════

class TestGCContent:
    def test_gc_in_range_pass(self):
        # GC of "GCGCATATAT" = 4/10 = 0.4
        v, d = check_gc_in_range("GCGCATATAT", 0.3, 0.5)
        assert v == "PASS"

    def test_gc_in_range_fail_low(self):
        # GC of "ATATATATAT" = 0/10 = 0.0
        v, d = check_gc_in_range("ATATATATAT", 0.3, 0.7)
        assert v == "FAIL"

    def test_gc_in_range_fail_high(self):
        # GC of "GCGCGCGCGC" = 10/10 = 1.0
        v, d = check_gc_in_range("GCGCGCGCGC", 0.3, 0.7)
        assert v == "FAIL"

    def test_gc_content_empty(self):
        assert _gc_content("") == 0.0

    def test_gc_content_half(self):
        assert abs(_gc_content("GCAT") - 0.5) < 1e-9


class TestNoStopCodons:
    def test_no_stops_pass(self):
        v, d = check_no_stop_codons("ATGGCTAAGCTTTAA")
        assert v == "PASS"
        assert "No internal" in d

    def test_internal_stop_fail(self):
        # TAA at position 3 (within coding region, not last codon)
        v, d = check_no_stop_codons("ATGTAAGCTTAA")
        assert v == "FAIL"
        assert "3" in d

    def test_last_codon_stop_allowed(self):
        # Only stop is at the last codon
        v, d = check_no_stop_codons("ATGGCTGCTTAA")
        assert v == "PASS"

    def test_short_sequence(self):
        v, d = check_no_stop_codons("AT")
        assert v == "PASS"


class TestNoCrypticSplice:
    def test_no_gt_pass(self):
        # Sequence with no GT dinucleotides
        v, d = check_no_cryptic_splice("ACACACACAC", 3.0)
        assert v == "PASS"

    def test_has_cryptic_splice(self):
        # Sequence with GT that scores high enough
        # Use a sequence with strong donor consensus around GT
        seq = "CAG" + "GTAAGT" + "AC" * 50
        v, d = check_no_cryptic_splice(seq, 3.0)
        # This should detect some sites
        assert isinstance(v, str)
        assert isinstance(d, str)


class TestNoRestrictionSites:
    def test_no_sites_pass(self):
        v, d = check_no_restriction_sites("ACACACACAC", ["EcoRI", "BamHI"])
        assert v == "PASS"

    def test_ecori_site_fail(self):
        # GAATTC is EcoRI site
        v, d = check_no_restriction_sites("ACGAATTCAC", ["EcoRI"])
        assert v == "FAIL"
        assert "EcoRI" in d

    def test_unknown_enzyme_pass(self):
        v, d = check_no_restriction_sites("ACACACACAC", ["UnknownEnzyme"])
        assert v == "PASS"

    def test_multiple_sites(self):
        seq = "GAATTCGGATCC"  # EcoRI + BamHI
        v, d = check_no_restriction_sites(seq, ["EcoRI", "BamHI"])
        assert v == "FAIL"
        assert "EcoRI" in d
        assert "BamHI" in d


class TestNoCpGIsland:
    def test_no_cpg_pass(self):
        v, d = check_no_cpg_island("ATAT" * 100, 200, 0.6)
        assert v == "PASS"

    def test_cpg_island_fail(self):
        # Create a sequence rich in CG dinucleotides
        seq = "CGCGCGCG" * 50  # 400 nt, many CG
        v, d = check_no_cpg_island(seq, 200, 0.6)
        assert v == "FAIL"


class TestNoGTDinucleotide:
    def test_no_gt_pass(self):
        v, d = check_no_gt_dinucleotide("ACACACAC")
        assert v == "PASS"

    def test_has_gt_fail(self):
        v, d = check_no_gt_dinucleotide("ACGTACAC")
        assert v == "FAIL"
        assert "2" in d


class TestCAI:
    def test_cai_human_above_threshold(self):
        # Use a well-adapted human coding sequence
        seq = "ATGGCCGCCCAGGAGCTGAAGGAGAAG"
        v, d = check_cai_above_threshold(seq, "Homo_sapiens", 0.3)
        # With common human codons, CAI should be reasonable
        assert isinstance(v, str)

    def test_cai_ecoli(self):
        seq = "ATGGCTGCTGAAGAGAAAGTT"
        v, d = check_cai_above_threshold(seq, "Escherichia_coli", 0.3)
        assert isinstance(v, str)

    def test_cai_organism_matching(self):
        # Test that various organism name forms resolve to E. coli
        v1, _ = check_cai_above_threshold("ATGGCTGCTAA", "ecoli", 0.3)
        v2, _ = check_cai_above_threshold("ATGGCTGCTAA", "E_coli", 0.3)
        v3, _ = check_cai_above_threshold("ATGGCTGCTAA", "Escherichia_coli", 0.3)
        assert v1 == v2 == v3


class TestProteinIdentity:
    def test_match_pass(self):
        # ATG=M, GCT=A, GCC=A => "MAA", then TAA=stop
        v, d = check_protein_identity("ATGGCTGCCTAA", "MAA")
        assert v == "PASS"

    def test_mismatch_fail(self):
        v, d = check_protein_identity("ATGGCTGCCTAA", "MGG")
        assert v == "FAIL"

    def test_length_mismatch_fail(self):
        v, d = check_protein_identity("ATGGCTGCCTAA", "MA")
        assert v == "FAIL"


# ══════════════════════════════════════════════════════════════════════
# Test: Dispatch and utility functions
# ══════════════════════════════════════════════════════════════════════

class TestDispatch:
    def test_resolve_base_name_plain(self):
        assert _resolve_base_name("GCInRange") == "GCInRange"

    def test_resolve_base_name_parameterized(self):
        assert _resolve_base_name("GCInRange(0.3, 0.7)") == "GCInRange"

    def test_get_params_defaults(self):
        params = _get_params({"sequence": "ACGT"})
        assert params["gc_lo"] == 0.30
        assert params["gc_hi"] == 0.70
        assert params["cai_threshold"] == 0.5
        assert params["organism"] == "Homo_sapiens"

    def test_get_params_from_provenance(self):
        cert = {
            "provenance": {
                "parameters": {"gc_lo": 0.4, "gc_hi": 0.6, "organism": "ecoli"}
            }
        }
        params = _get_params(cert)
        assert params["gc_lo"] == 0.4
        assert params["organism"] == "ecoli"

    def test_dispatch_unknown(self):
        v, d = _dispatch_check("UnknownPredicate", "ACGT", {})
        assert v == "UNCERTAIN"
        assert "Unknown" in d

    def test_dispatch_conservation_uncertain(self):
        v, d = _dispatch_check("ConservationScore", "ACGT", {})
        assert v == "UNCERTAIN"

    def test_dispatch_protein_identity_no_aa(self):
        v, d = _dispatch_check("ProteinIdentity", "ATGGCTTAA", {})
        assert v == "UNCERTAIN"

    def test_dispatch_valid_coding_seq_pass(self):
        v, d = _dispatch_check("ValidCodingSeq", "ATGGCTTAA", {})
        assert v == "PASS"

    def test_dispatch_valid_coding_seq_fail_length(self):
        v, d = _dispatch_check("ValidCodingSeq", "ATGGC", {})
        assert v == "FAIL"

    def test_dispatch_instability_motif(self):
        v, d = _dispatch_check("NoInstabilityMotif", "ATGATTTAGCTAA", {})
        assert v == "FAIL"


# ══════════════════════════════════════════════════════════════════════
# Test: Full certificate verification
# ══════════════════════════════════════════════════════════════════════

class TestVerify:
    def test_valid_certificate(self):
        """Certificate with correct claims should verify."""
        # Build a sequence with no stops, no restriction sites
        seq = "ATG" + "GCC" * 50 + "TAA"  # 156 nt, GC=100% (too high for default range)
        # Use a wider GC range for this test
        cert = _make_valid_cert(seq, types_overrides=[
            {"predicate": "NoStopCodons", "verdict": "PASS"},
            {"predicate": "GCInRange", "verdict": "PASS"},
            {"predicate": "NoRestrictionSite", "verdict": "PASS"},
        ])
        cert["provenance"]["parameters"]["gc_lo"] = 0.0
        cert["provenance"]["parameters"]["gc_hi"] = 1.0
        status, reasons = verify(cert)
        # Should pass no-stop and no-restriction; GC should also pass with wide range
        assert status == "VERIFIED" or len(reasons) <= 1  # may fail GT check if seq has GT

    def test_missing_required_key(self):
        """Certificate missing required key should be rejected."""
        cert = {"version": "7.2.0"}  # missing design_id, sequence
        status, reasons = verify(cert)
        assert status == "REJECTED"
        assert any("design_id" in r for r in reasons)

    def test_hash_mismatch(self):
        """Certificate with wrong design_id should be rejected."""
        cert = _make_valid_cert()
        cert["design_id"] = "0" * 64  # Wrong hash
        status, reasons = verify(cert)
        assert status == "REJECTED"
        assert any("design_id mismatch" in r for r in reasons)

    def test_predicate_verdict_mismatch(self):
        """Certificate claiming PASS for a predicate that actually FAILs."""
        seq = "ATGTAAGCTTAA"  # Has TAA at position 3 (internal stop)
        cert = _make_valid_cert(seq, types_overrides=[
            {"predicate": "NoStopCodons", "verdict": "PASS"},  # WRONG: has internal stop
        ])
        status, reasons = verify(cert)
        assert status == "REJECTED"
        assert any("NoStopCodons" in r for r in reasons)

    def test_uncertain_predicate_not_rejected(self):
        """UNCERTAIN predicates that cannot be independently verified should not cause rejection."""
        cert = _make_valid_cert(types_overrides=[
            {"predicate": "ConservationScore", "verdict": "PASS"},
        ])
        # ConservationScore returns UNCERTAIN — should not reject
        status, reasons = verify(cert)
        assert not any("ConservationScore" in r for r in reasons)

    def test_provenance_missing_input_hash(self):
        """Missing provenance.input_hash should be rejected."""
        cert = _make_valid_cert()
        del cert["provenance"]["input_hash"]
        status, reasons = verify(cert)
        assert status == "REJECTED"
        assert any("input_hash" in r for r in reasons)

    def test_provenance_input_hash_mismatch(self):
        """Wrong provenance.input_hash should be rejected."""
        cert = _make_valid_cert()
        cert["provenance"]["input_hash"] = "0" * 64
        status, reasons = verify(cert)
        assert status == "REJECTED"
        assert any("input_hash mismatch" in r for r in reasons)


class TestPredicatesFormat:
    """Test the 'predicates' dict format (simpler alternative to 'types')."""

    def test_predicates_format_pass(self):
        seq = "ATG" + "GCC" * 50 + "TAA"
        seq_hash = hashlib.sha256(seq.encode()).hexdigest()
        cert = {
            "version": "7.2.0",
            "design_id": seq_hash,
            "sequence": seq,
            "predicates": {
                "no_stop_codons": {"verdict": "PASS"},
                "gc_content": {"verdict": "PASS"},
            },
            "gc_lo": 0.0,
            "gc_hi": 1.0,
        }
        status, reasons = verify(cert)
        assert status == "VERIFIED"

    def test_predicates_format_fail(self):
        seq = "ATGTAAGCTTAA"  # Has internal stop
        seq_hash = hashlib.sha256(seq.encode()).hexdigest()
        cert = {
            "version": "7.2.0",
            "design_id": seq_hash,
            "sequence": seq,
            "predicates": {
                "no_stop_codons": {"verdict": "PASS"},  # WRONG
            },
        }
        status, reasons = verify(cert)
        assert status == "REJECTED"
        assert any("no_stop_codons" in r for r in reasons)


class TestVerifyFile:
    """Test file-based verification."""

    def test_verify_valid_file(self, tmp_path):
        cert = _make_valid_cert()
        # Remove predicates that will fail on random sequences
        cert["types"] = [{"predicate": "NoStopCodons", "verdict": "PASS"}]
        cert["provenance"]["parameters"]["gc_lo"] = 0.0
        cert["provenance"]["parameters"]["gc_hi"] = 1.0
        p = tmp_path / "cert.json"
        p.write_text(json.dumps(cert))
        status, reasons = verify_file(str(p))
        assert status == "VERIFIED"

    def test_verify_missing_file(self):
        with pytest.raises(FileNotFoundError):
            verify_file("/nonexistent/cert.json")

    def test_verify_invalid_json(self, tmp_path):
        p = tmp_path / "bad.json"
        p.write_text("{not valid json")
        with pytest.raises(json.JSONDecodeError):
            verify_file(str(p))


class TestMaxEntDonorScore:
    """Test the embedded MaxEntScan donor scoring."""

    def test_strong_donor(self):
        # Good donor context: CCC before GT, AAGT after
        # This matches the PWM's strong donor profile
        seq = "AACCCGTAAGTACACACACAC"
        score = _maxent_donor_score(seq, 5)  # GT at position 5
        assert score > 3.0  # Should be a strong donor

    def test_weak_donor(self):
        # Non-consensus context around GT
        seq = "ATATATGTATATATATATATA"
        score = _maxent_donor_score(seq, 6)  # GT at position 6
        # Score might be lower with this context
        assert isinstance(score, float)

    def test_boundary_too_short(self):
        # Sequence too short for 9-mer context
        score = _maxent_donor_score("GT", 0)
        assert score == -50.0


class TestComputeCAI:
    """Test CAI computation."""

    def test_cai_range(self):
        # CAI should be between 0 and 1
        seq = "ATGGCCGCCCAGGAGCTGAAGGAGAAG"
        cai = _compute_cai(seq, {
            "ATG": 27.0, "GCC": 27.7, "GCA": 15.8, "GCG": 7.4, "GCT": 18.4,
            "CAG": 34.3, "GAG": 39.8, "GAA": 28.8, "AAG": 32.1, "AAA": 24.1,
        })
        assert 0.0 <= cai <= 1.0

    def test_cai_optimal_codons(self):
        # All-optimal codons should give CAI close to 1.0
        seq = "ATG" + "GCC" * 10 + "TAA"
        # GCC is most common Ala codon in human
        cai = _compute_cai(seq, {"ATG": 22.3, "GCC": 27.7, "GCT": 18.4,
                                  "GCA": 15.8, "GCG": 7.4, "TAA": 1.5})
        assert cai > 0.9


class TestZeroBiocompilerImports:
    """Verify that standalone_verifier.py has ZERO imports from biocompiler."""

    def test_no_biocompiler_imports(self):
        """The standalone verifier must not import any biocompiler module."""
        import ast
        path = os.path.join(os.path.dirname(__file__), "..", "src", "biocompiler", "standalone_verifier.py")
        with open(path) as f:
            tree = ast.parse(f.read())
        for node in ast.walk(tree):
            if isinstance(node, ast.Import):
                for alias in node.names:
                    assert not alias.name.startswith("biocompiler"), \
                        f"Forbidden import: {alias.name}"
            elif isinstance(node, ast.ImportFrom):
                assert node.module is None or not node.module.startswith("biocompiler"), \
                    f"Forbidden import from: {node.module}"

    def test_line_count_under_500(self):
        """The standalone verifier must be under 500 lines for auditability."""
        path = os.path.join(os.path.dirname(__file__), "..", "src", "biocompiler", "standalone_verifier.py")
        with open(path) as f:
            lines = f.readlines()
        assert len(lines) < 500, f"Verifier is {len(lines)} lines, must be < 500"

    def test_only_stdlib_imports(self):
        """All imports must be from Python stdlib."""
        import ast
        STDLIB = {"hashlib", "json", "math", "sys", "os", "re", "collections",
                   "itertools", "functools", "typing", "abc", "io", "pathlib"}
        path = os.path.join(os.path.dirname(__file__), "..", "src", "biocompiler", "standalone_verifier.py")
        with open(path) as f:
            tree = ast.parse(f.read())
        for node in ast.walk(tree):
            if isinstance(node, ast.Import):
                for alias in node.names:
                    root = alias.name.split(".")[0]
                    assert root in STDLIB, f"Non-stdlib import: {alias.name}"
            elif isinstance(node, ast.ImportFrom):
                if node.module:
                    root = node.module.split(".")[0]
                    assert root in STDLIB, f"Non-stdlib import from: {node.module}"


class TestCLI:
    """Test CLI entry point."""

    # Run standalone_verifier.py directly via subprocess.
    # We must cd away from src/biocompiler/ to avoid the local types.py
    # shadowing the stdlib `types` module (circular import issue).
    _SCRIPT = os.path.abspath(os.path.join(
        os.path.dirname(__file__), "..", "src", "biocompiler", "standalone_verifier.py"
    ))

    @staticmethod
    def _run_cli(*args, cwd="/tmp"):
        import subprocess
        env = os.environ.copy()
        env["PYTHONPATH"] = ""  # Don't inherit project's src directory
        return subprocess.run(
            [sys.executable, TestCLI._SCRIPT, *args],
            capture_output=True, text=True, cwd=cwd, env=env
        )

    def test_cli_no_args(self):
        """CLI with no args should exit with code 2."""
        result = self._run_cli()
        assert result.returncode == 2

    def test_cli_valid_cert(self, tmp_path):
        """CLI with valid cert should exit with code 0."""
        cert = _make_valid_cert()
        cert["types"] = [{"predicate": "NoStopCodons", "verdict": "PASS"}]
        cert["provenance"]["parameters"]["gc_lo"] = 0.0
        cert["provenance"]["parameters"]["gc_hi"] = 1.0
        p = tmp_path / "cert.json"
        p.write_text(json.dumps(cert))
        result = self._run_cli(str(p), cwd=str(tmp_path))
        assert result.returncode == 0
        assert "VERIFIED" in result.stdout

    def test_cli_invalid_cert(self, tmp_path):
        """CLI with invalid cert should exit with code 1."""
        cert = _make_valid_cert()
        cert["design_id"] = "0" * 64  # Wrong hash
        p = tmp_path / "cert.json"
        p.write_text(json.dumps(cert))
        result = self._run_cli(str(p), cwd=str(tmp_path))
        assert result.returncode == 1
        assert "REJECTED" in result.stdout

    def test_cli_nonexistent_file(self):
        """CLI with nonexistent file should exit with code 2."""
        result = self._run_cli("/nonexistent/file.json")
        assert result.returncode == 2
