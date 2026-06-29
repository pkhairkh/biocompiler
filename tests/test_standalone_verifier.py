"""TIGHTEN-5 — End-to-end tests for the standalone certificate verifier.

The paper claims BioCompiler ships a "standalone verifier (~460 LOC)" that
re-checks every predicate from a serialized certificate independently of
the BioCompiler codebase (stdlib-only).  These tests pin that claim by:

  1. Verifying the script exists and is roughly the claimed LOC count.
  2. Round-tripping the verifier's own ``--generate-test`` certificate.
  3. Building a verifier-format certificate from a real optimizer run
     and confirming the verifier re-checks all 14 core predicates.
  4. Confirming the verifier detects a tampered design_id (hash integrity).
  5. Confirming the verifier rejects a certificate with a missing required
     field.
  6. Confirming the verifier can read a certificate from stdin.

The verifier is invoked as a subprocess so that the test exercises the
real CLI entry point (``python scripts/standalone_verifier.py <cert>``)
exactly as an external consumer would.
"""
from __future__ import annotations

import json
import os
import subprocess
import sys
import tempfile

import pytest

# The standalone verifier lives outside the importable package.
# Resolve relative to the test file so the test runs from any checkout.
_VERIFY_SCRIPT = os.path.abspath(
    os.path.join(os.path.dirname(__file__), "..", "scripts", "standalone_verifier.py")
)
_PYTHON = sys.executable
_SCRIPTS_DIR = os.path.dirname(_VERIFY_SCRIPT)


def _run_verifier(*args: str, timeout: int = 30) -> subprocess.CompletedProcess:
    """Invoke the standalone verifier as a subprocess."""
    return subprocess.run(
        [_PYTHON, _VERIFY_SCRIPT, *args],
        capture_output=True,
        text=True,
        timeout=timeout,
    )


def _write_temp_cert(cert: dict) -> str:
    """Write a certificate dict to a temp JSON file; return its path."""
    fd, path = tempfile.mkstemp(suffix=".json", prefix="biocompiler_cert_")
    with os.fdopen(fd, "w") as f:
        json.dump(cert, f)
    return path


def _import_standalone_verifier():
    """Import the standalone_verifier module in-process for helper checks."""
    if _SCRIPTS_DIR not in sys.path:
        sys.path.insert(0, _SCRIPTS_DIR)
    import standalone_verifier as sv  # type: ignore[import-not-found]
    return sv


def _build_standalone_cert_from_optimizer(
    protein: str = "MVHLTPEEK",
    organism: str = "human",
) -> dict:
    """Run the optimizer and translate its output into the standalone
    verifier's certificate JSON format.

    The standalone verifier expects keys like ``dna_sequence`` and
    ``protein_sequence`` (see the module docstring of
    ``scripts/standalone_verifier.py``).  The in-package Certificate
    dataclass uses ``sequence`` / ``types`` instead, so we adapt.

    The design_id is recomputed with the verifier's own hash function so
    that the verifier, on read, will compute the exact same hash.
    """
    sv = _import_standalone_verifier()
    from biocompiler.optimizer.pipeline_core import optimize_sequence
    from biocompiler.provenance.certificate import compute_certificate

    result = optimize_sequence(protein, organism=organism, strict_mode=False)
    level = compute_certificate(result.predicate_results)

    org_field = "Homo_sapiens" if organism in ("human", "Homo_sapiens") else organism
    cert = {
        "design_id": "",
        "organism": org_field,
        "protein_sequence": result.protein,
        "dna_sequence": result.sequence,
        "original_protein": result.protein,
        "original_dna_sequence": result.sequence,
        "gc_range": [0.30, 0.70],
        "certificate_level": level.value,  # GOLD / SILVER / BRONZE
        "predicates": {},
        "optimizer_version": "1.0.0",
    }
    cert["design_id"] = sv.compute_design_id(cert)
    return cert


class TestStandaloneVerifierExists:
    """Structural checks: the verifier exists and is ~460 LOC."""

    def test_script_exists(self):
        assert os.path.exists(_VERIFY_SCRIPT), (
            f"standalone_verifier.py not found at {_VERIFY_SCRIPT}"
        )

    def test_script_loc_matches_paper_claim(self):
        """The paper advertises ~460 LOC; allow a small window around it."""
        with open(_VERIFY_SCRIPT) as f:
            loc = sum(1 for _ in f)
        # 460 claimed; tolerate +/-10% drift for header/footer churn.
        assert 410 <= loc <= 520, f"Standalone verifier LOC={loc} outside expected ~460"

    def test_script_has_cli_main(self):
        """The script must define a ``main`` entry point and a
        ``verify_certificate`` function so it can be used both as a CLI
        and as an importable module."""
        with open(_VERIFY_SCRIPT) as f:
            src = f.read()
        assert "def main(" in src
        assert "def verify_certificate(" in src
        assert 'if __name__ == "__main__"' in src

    def test_script_uses_only_stdlib(self):
        """The script must rely only on the Python standard library so
        that it can be shipped standalone (no BioCompiler dependency)."""
        with open(_VERIFY_SCRIPT) as f:
            src = f.read()
        import_lines = [
            ln.strip()
            for ln in src.splitlines()
            if ln.strip().startswith(("import ", "from "))
        ]
        assert import_lines, "No import statements found"
        forbidden_prefixes = ("biocompiler", "numpy", "scipy", "Bio", "pandas")
        for ln in import_lines:
            # Extract the module name being imported.
            tokens = ln.split()
            module = tokens[1] if tokens[0] == "import" else tokens[1]
            # Strip any "as"/"," suffixes.
            module = module.split(",")[0].split(" as ")[0]
            for bad in forbidden_prefixes:
                assert not module.startswith(bad), (
                    f"Standalone verifier imports non-stdlib dependency: {ln}"
                )


class TestStandaloneVerifierGenerateTestCertificate:
    """``--generate-test`` produces a syntactically valid certificate."""

    def test_generate_test_outputs_valid_json(self):
        proc = _run_verifier("--generate-test")
        assert proc.returncode == 0, proc.stderr
        cert = json.loads(proc.stdout)
        for key in (
            "design_id",
            "organism",
            "protein_sequence",
            "dna_sequence",
            "gc_range",
            "certificate_level",
        ):
            assert key in cert, f"--generate-test output missing key: {key}"

    def test_generate_test_certificate_level_is_gold(self):
        """The bundled test certificate is advertised as GOLD."""
        proc = _run_verifier("--generate-test")
        cert = json.loads(proc.stdout)
        assert cert["certificate_level"] == "GOLD"

    def test_generate_test_design_id_is_self_consistent(self):
        """The design_id written by --generate-test must match a fresh
        recompute using the verifier's own hash function."""
        sv = _import_standalone_verifier()
        proc = _run_verifier("--generate-test")
        cert = json.loads(proc.stdout)
        assert cert["design_id"] == sv.compute_design_id(cert)


class TestStandaloneVerifierEndToEnd:
    """End-to-end: optimizer output -> certificate JSON -> standalone verifier."""

    def test_verifier_accepts_optimizer_generated_certificate(self):
        """A certificate built from a real optimizer run must verify
        cleanly: every one of the 14 core predicates must PASS and the
        overall verdict must be PASS (exit code 0).

        This pins the paper's claim that the optimizer produces
        sequences the standalone verifier can independently confirm.
        """
        cert = _build_standalone_cert_from_optimizer()
        cert_path = _write_temp_cert(cert)
        try:
            proc = _run_verifier(cert_path)
        finally:
            os.unlink(cert_path)

        assert proc.returncode == 0, (
            f"Verifier rejected optimizer certificate.\n"
            f"stdout:\n{proc.stdout}\nstderr:\n{proc.stderr}"
        )
        assert "Core predicates: 14 PASS / 0 FAIL / 14 total" in proc.stdout
        assert "Overall verdict : PASS" in proc.stdout
        assert "_hash_integrity" in proc.stdout
        assert "SLOT predicates: 22 UNCERTAIN" in proc.stdout

    def test_verifier_independently_rechecks_predicates(self):
        """The verifier must NOT trust the certificate's claimed level
        blindly -- it re-runs every predicate from first principles.

        We verify this by giving it a certificate that *claims* GOLD but
        whose DNA sequence actually contains a CpG dinucleotide (which
        the NoCpGIsland predicate flags as FAIL).  The verifier should
        detect this and exit non-zero.
        """
        sv = _import_standalone_verifier()

        # ATG GTG CAC CGC TGA -> Met-Val-His-Arg-stop.  Contains "CG"
        # at position 9-10 (CpG) which NoCpGIsland flags as FAIL.
        bad_dna = "ATGGTGCACCGCTGA"
        cert = {
            "design_id": "",
            "organism": "Homo_sapiens",
            "protein_sequence": "MVHR",
            "dna_sequence": bad_dna,
            "original_protein": "MVHR",
            "original_dna_sequence": bad_dna,
            "gc_range": [0.30, 0.70],
            "certificate_level": "GOLD",  # dishonest claim
            "predicates": {},
            "optimizer_version": "1.0.0",
        }
        cert["design_id"] = sv.compute_design_id(cert)

        cert_path = _write_temp_cert(cert)
        try:
            proc = _run_verifier(cert_path)
        finally:
            os.unlink(cert_path)

        assert proc.returncode == 1, (
            f"Verifier should have rejected the dishonest GOLD certificate.\n"
            f"stdout:\n{proc.stdout}"
        )
        assert "NoCpGIsland" in proc.stdout
        assert "FAIL" in proc.stdout
        assert "Overall verdict : FAIL" in proc.stdout

    def test_verifier_detects_tampered_design_id(self):
        """If the design_id hash is mutated, the verifier must flag the
        ``_hash_integrity`` predicate as FAIL.

        Note: the standalone verifier reports hash integrity separately
        from the core predicate count.  The standalone CLI's overall
        PASS/FAIL verdict is driven only by the 14 core predicate
        verdicts (see ``print_results`` in ``standalone_verifier.py``),
        so a tampered hash surfaces as a FAIL line in the report rather
        than a non-zero exit code.  We assert both: the FAIL line is
        present, and the expected-vs-actual hashes are printed.
        """
        cert = _build_standalone_cert_from_optimizer()
        bogus_hash = "sha256:" + "0" * 64
        cert["design_id"] = bogus_hash
        cert_path = _write_temp_cert(cert)
        try:
            proc = _run_verifier(cert_path)
        finally:
            os.unlink(cert_path)

        # The verifier must surface the hash mismatch in its report.
        assert "_hash_integrity" in proc.stdout
        assert "FAIL" in proc.stdout
        # The bogus hash we wrote must appear, alongside the recomputed
        # expected hash.
        assert "0000000000000000000000000000000000000000000000000000000000000000" in proc.stdout
        assert "Expected sha256:" in proc.stdout

    def test_verifier_reads_from_stdin(self):
        """The verifier supports ``-`` as a stdin path (per its CLI doc)."""
        cert = _build_standalone_cert_from_optimizer()
        proc = subprocess.run(
            [_PYTHON, _VERIFY_SCRIPT, "-"],
            input=json.dumps(cert),
            capture_output=True,
            text=True,
            timeout=30,
        )
        assert proc.returncode == 0, (
            f"stdin verification failed.\nstdout:\n{proc.stdout}\nstderr:\n{proc.stderr}"
        )
        assert "Overall verdict : PASS" in proc.stdout

    def test_verifier_rejects_missing_required_field(self):
        """A certificate missing ``dna_sequence`` must be rejected with
        a non-zero exit code (the script documents exit code 3 for this
        case)."""
        cert = _build_standalone_cert_from_optimizer()
        del cert["dna_sequence"]
        cert_path = _write_temp_cert(cert)
        try:
            proc = _run_verifier(cert_path)
        finally:
            os.unlink(cert_path)

        assert proc.returncode != 0
        assert proc.returncode == 3
        assert "missing required fields" in proc.stderr.lower()

    def test_verifier_rejects_nonexistent_file(self):
        """A non-existent path must exit non-zero (documented exit code 2)."""
        proc = _run_verifier("/tmp/does_not_exist_biocompiler_cert.json")
        assert proc.returncode == 2
        assert "not found" in proc.stderr.lower()

    def test_verifier_no_args_prints_usage(self):
        """With no args, the verifier should print its docstring/usage."""
        proc = _run_verifier()
        assert proc.returncode == 1
        assert "Usage" in proc.stderr or "standalone_verifier" in proc.stderr


class TestStandaloneVerifierPredicateCoverage:
    """The verifier must re-check every predicate the paper claims.

    The script's docstring advertises:
        14 CORE predicates -- fully re-checked from first principles
        22 SLOT  predicates -- marked UNCERTAIN (external oracles required)
    """

    def test_has_14_core_predicate_checks(self):
        sv = _import_standalone_verifier()
        cert = _build_standalone_cert_from_optimizer()
        results = sv.verify_certificate(cert)
        core_names = {
            n for n in results
            if not n.startswith("_") and n not in sv.SLOT_PREDICATES
        }
        assert len(core_names) == 14, (
            f"Expected 14 core predicates, got {len(core_names)}: "
            f"{sorted(core_names)}"
        )

    def test_has_22_slot_predicates(self):
        sv = _import_standalone_verifier()
        cert = _build_standalone_cert_from_optimizer()
        results = sv.verify_certificate(cert)
        slot_results = {
            n: r for n, r in results.items()
            if n in sv.SLOT_PREDICATES
        }
        assert len(slot_results) == 22, (
            f"Expected 22 SLOT predicates, got {len(slot_results)}"
        )
        for name, res in slot_results.items():
            assert res["verdict"] == "UNCERTAIN", (
                f"SLOT predicate {name} should be UNCERTAIN, got {res['verdict']}"
            )

    def test_hash_integrity_predicate_always_present(self):
        sv = _import_standalone_verifier()
        cert = _build_standalone_cert_from_optimizer()
        results = sv.verify_certificate(cert)
        assert "_hash_integrity" in results
        assert results["_hash_integrity"]["verdict"] == "PASS"


if __name__ == "__main__":
    import pytest
    sys.exit(pytest.main([__file__, "-v"]))
