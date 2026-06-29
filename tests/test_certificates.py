"""TIGHTEN-5 — Graduated certificate tests (GOLD / SILVER / BRONZE).

The paper claims BioCompiler ships "graduated certificates (Gold/Silver/
Bronze)" that document what was verified, not just whether everything
passed.  These tests pin that claim by exercising:

  1. ``compute_certificate`` — the tier-computation rules.
  2. ``format_certificate`` — the human-readable report.
  3. ``generate_certificate`` — the full Certificate dataclass builder.
  4. ``verify_certificate`` — the in-package independent re-verifier.
  5. End-to-end: optimizer -> predicate_results -> certificate tier.

Tier rules (from ``certificate.compute_certificate``):
    GOLD   : all predicates PASS, 0 UNCERTAIN, no mutagenesis/unavoidable
    SILVER : 1 UNCERTAIN OR mutagenesis OR unavoidable constraints
    BRONZE : >=2 UNCERTAIN OR any predicate not satisfied

Uncertainty capping (Issue #10):
    * LIKELY_PASS / LIKELY_FAIL do NOT cap (calibrated uncertainty).
    * Only true UNCERTAIN verdicts cap the tier.
"""
from __future__ import annotations

import os
import sys

import pytest

# The standalone verifier lives outside the package; needed for the
# end-to-end "standalone cert from optimizer" test. Resolve relative to
# this test file so the test runs from any checkout.
_SCRIPTS_DIR = os.path.abspath(
    os.path.join(os.path.dirname(__file__), "..", "scripts")
)
if _SCRIPTS_DIR not in sys.path:
    sys.path.insert(0, _SCRIPTS_DIR)

from biocompiler.optimizer.pipeline_core import optimize_sequence
from biocompiler.provenance.certificate import (
    compute_certificate,
    compute_uncertainty_summary,
    format_certificate,
    generate_certificate,
    verify_certificate,
)
from biocompiler.shared.types import SLOTMode, TypeCheckResult, Verdict
from biocompiler.type_system import CertLevel, PredicateResult, SpliceVerdict


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def _all_pass_results() -> list[PredicateResult]:
    """A clean predicate-result list where every predicate PASSes with
    no UNCERTAIN, no mutagenesis, no unavoidable constraints."""
    return [
        PredicateResult("NoStopCodons", True, verdict=Verdict.PASS,
                        details="No internal stop codons"),
        PredicateResult("NoCrypticSplice", True, verdict=SpliceVerdict.PASS,
                        details="No GT dinucleotides found"),
        PredicateResult("NoCpGIsland", True, verdict=Verdict.PASS,
                        details="Worst CpG Obs/Exp ratio 0.300 <= 0.6"),
        PredicateResult("NoRestrictionSite", True, verdict=Verdict.PASS,
                        details="No restriction sites found"),
        PredicateResult("NoGTDinucleotide", True, verdict=Verdict.PASS,
                        details="No GT dinucleotides found"),
        PredicateResult("ValidCodingSeq", True, verdict=Verdict.PASS,
                        details="All codons valid"),
        PredicateResult("ConservationScore", True, verdict=Verdict.PASS,
                        details="All AA conservation scores >= -1"),
        PredicateResult("CodonOptimality", True, verdict=Verdict.PASS,
                        details="Worst CAI: GCT=0.7244, min=0.0"),
    ]


# ---------------------------------------------------------------------------
# CertLevel enum
# ---------------------------------------------------------------------------

class TestCertLevelEnum:
    """The paper claims exactly three tiers: Gold, Silver, Bronze."""

    def test_certlevel_has_exactly_three_tiers(self):
        tiers = {level.value for level in CertLevel}
        assert tiers == {"GOLD", "SILVER", "BRONZE"}, (
            f"Expected exactly GOLD/SILVER/BRONZE, got {sorted(tiers)}"
        )

    def test_certlevel_ordering_is_gold_gt_silver_gt_bronze(self):
        """GOLD is the strongest tier; BRONZE the weakest.

        We rely on Python's value-based enum comparison by name to make
        sure nobody accidentally reorders the tiers.
        """
        names = [level.name for level in CertLevel]
        assert names.index("GOLD") < names.index("SILVER") < names.index("BRONZE")


# ---------------------------------------------------------------------------
# compute_certificate tier rules
# ---------------------------------------------------------------------------

class TestComputeCertificateTiers:
    """Pin the GOLD / SILVER / BRONZE decision rules."""

    def test_gold_certificate(self):
        """All predicates satisfied with no mutagenesis -> GOLD."""
        results = _all_pass_results()
        cert = compute_certificate(results)
        assert cert == CertLevel.GOLD
        assert cert.value == "GOLD"

    def test_silver_certificate_unavoidable(self):
        """All predicates passed but some have unavoidable constraints
        -> SILVER."""
        results = _all_pass_results()
        # Mark one predicate as having unavoidable constraints.
        results[4] = PredicateResult(
            "NoGTDinucleotide", True, verdict=Verdict.PASS,
            details="All 2 GT dinucleotides are unavoidable",
        )
        cert = compute_certificate(results)
        assert cert == CertLevel.SILVER

    def test_silver_certificate_mutagenesis(self):
        """All predicates passed but some involved mutagenesis -> SILVER."""
        results = _all_pass_results()
        results[4] = PredicateResult(
            "NoGTDinucleotide", True, verdict=Verdict.PASS,
            details="No GT dinucleotides found mutagenesis applied: pos 3:V->I",
        )
        cert = compute_certificate(results)
        assert cert == CertLevel.SILVER

    def test_silver_certificate_structured_mutagenesis_flag(self):
        """The structured ``mutagenesis_applied=True`` flag also yields
        SILVER (preferred over string-matching per the docstring)."""
        results = _all_pass_results()
        results[0] = PredicateResult(
            "NoStopCodons", True, verdict=Verdict.PASS, details="ok",
            mutagenesis_applied=True,
        )
        cert = compute_certificate(results)
        assert cert == CertLevel.SILVER

    def test_bronze_certificate_unsatisfied(self):
        """Any unsatisfied predicate -> BRONZE."""
        results = _all_pass_results()
        results[4] = PredicateResult(
            "NoGTDinucleotide", False, verdict=Verdict.FAIL,
            details="Avoidable GT dinucleotides at [3]",
        )
        cert = compute_certificate(results)
        assert cert == CertLevel.BRONZE

    def test_silver_certificate_one_uncertain(self):
        """A single UNCERTAIN verdict caps the tier at SILVER
        (uncertainty capping rule, Issue #10)."""
        results = _all_pass_results()
        results.append(PredicateResult(
            "mRNASecondaryStructure", True, verdict=Verdict.UNCERTAIN,
            details="Nussinov fallback: dG~0.0 kcal/mol (uncertain)",
        ))
        cert = compute_certificate(results)
        assert cert == CertLevel.SILVER

    def test_bronze_certificate_two_uncertain(self):
        """Two or more UNCERTAIN verdicts cap the tier at BRONZE."""
        results = _all_pass_results()
        results.append(PredicateResult(
            "mRNASecondaryStructure", True, verdict=Verdict.UNCERTAIN,
            details="uncertain",
        ))
        results.append(PredicateResult(
            "ProteinSolubility", True, verdict=Verdict.UNCERTAIN,
            details="uncertain",
        ))
        cert = compute_certificate(results)
        assert cert == CertLevel.BRONZE

    def test_likely_pass_does_not_cap_to_silver(self):
        """LIKELY_PASS expresses *calibrated* uncertainty and must NOT
        cap the tier (per Issue #10: only true UNCERTAIN caps)."""
        results = _all_pass_results()
        results.append(PredicateResult(
            "SomeHeuristic", True, verdict=Verdict.LIKELY_PASS,
            details="likely ok",
        ))
        cert = compute_certificate(results)
        assert cert == CertLevel.GOLD, (
            f"LIKELY_PASS must not cap; got {cert}"
        )

    def test_likely_fail_does_not_cap(self):
        """LIKELY_FAIL is treated as a calibrated verdict, not as
        UNCERTAIN; it does not cap to SILVER by itself.

        Note: ``passed=False`` for LIKELY_FAIL still triggers BRONZE via
        the unsatisfied-predicate rule.  Here we set ``passed=True`` to
        isolate the uncertainty-capping rule.
        """
        results = _all_pass_results()
        results.append(PredicateResult(
            "SomeHeuristic", True, verdict=Verdict.LIKELY_FAIL,
            details="likely fail but treated as passing for test",
        ))
        cert = compute_certificate(results)
        # Not UNCERTAIN, so does not cap to SILVER.  Stays GOLD.
        assert cert == CertLevel.GOLD


# ---------------------------------------------------------------------------
# compute_uncertainty_summary
# ---------------------------------------------------------------------------

class TestUncertaintySummary:
    """``compute_uncertainty_summary`` reports the verdict distribution."""

    def test_summary_counts_verdicts(self):
        """Build a result list with a known verdict distribution and
        confirm the summary counts match.

        We use plain ``Verdict`` values (not ``SpliceVerdict``) so that
        every verdict matches one of the summary's branches.
        """
        results = [
            PredicateResult("A", True, verdict=Verdict.PASS, details=""),
            PredicateResult("B", True, verdict=Verdict.PASS, details=""),
            PredicateResult("C", True, verdict=Verdict.UNCERTAIN, details="?"),
            PredicateResult("D", True, verdict=Verdict.LIKELY_PASS, details="?"),
        ]
        summary = compute_uncertainty_summary(results)
        assert summary["total_predicates"] == 4
        assert summary["uncertain_count"] == 1
        assert "C" in summary["uncertain_predicates"]
        assert summary["likely_pass_count"] == 1
        assert summary["likely_fail_count"] == 0
        assert summary["definite_count"] == 2

    def test_summary_confidence_score_in_range(self):
        """Confidence is a weighted average in [0.0, 1.0]."""
        results = [
            PredicateResult("A", True, verdict=Verdict.PASS, details=""),
            PredicateResult("B", True, verdict=Verdict.PASS, details=""),
        ]
        summary = compute_uncertainty_summary(results)
        # All PASS -> confidence 1.0.
        assert summary["confidence_score"] == 1.0
        assert 0.0 <= summary["confidence_score"] <= 1.0

    def test_summary_confidence_score_with_uncertain(self):
        """A mix of PASS + UNCERTAIN yields a fractional confidence."""
        results = [
            PredicateResult("A", True, verdict=Verdict.PASS, details=""),
            PredicateResult("B", True, verdict=Verdict.UNCERTAIN, details="?"),
        ]
        summary = compute_uncertainty_summary(results)
        # (1.0 + 0.5) / 2 = 0.75
        assert summary["confidence_score"] == 0.75


# ---------------------------------------------------------------------------
# format_certificate
# ---------------------------------------------------------------------------

class TestFormatCertificate:
    """``format_certificate`` produces a human-readable report."""

    def test_format_certificate_contains_level(self):
        results = _all_pass_results()
        text = format_certificate(results, "ATGGCTTAA", "ecoli")
        assert "GOLD" in text

    def test_format_certificate_contains_predicate(self):
        results = [
            PredicateResult("NoStopCodons", True, verdict=Verdict.PASS,
                            details="No internal stop codons"),
            PredicateResult("ValidCodingSeq", False, verdict=Verdict.FAIL,
                            details="Sequence length 8 not divisible by 3"),
        ]
        text = format_certificate(results, "ATGGCTGC", "ecoli")
        assert "NoStopCodons" in text
        assert "ValidCodingSeq" in text
        assert "PASS" in text
        assert "FAIL" in text

    def test_format_certificate_describes_all_three_tiers(self):
        """The report footer must explain what each tier means so that
        downstream consumers can interpret the level."""
        results = _all_pass_results()
        text = format_certificate(results, "ATGGCTTAA", "ecoli")
        # The footer lines describe all three tiers.
        assert "GOLD" in text
        assert "SILVER" in text
        assert "BRONZE" in text

    def test_format_certificate_marks_slot_predicates(self):
        """SLOT predicates must be marked with an asterisk in the report
        so consumers can distinguish formal vs empirical evidence."""
        # NoCrypticSplice is a SLOT predicate in the BioCompiler registry.
        results = _all_pass_results()
        text = format_certificate(results, "ATGGCTTAA", "ecoli")
        assert "SLOT" in text

    def test_format_certificate_includes_version(self):
        """The report header must show the BioCompiler version."""
        from biocompiler.provenance.certificate import VERSION
        results = _all_pass_results()
        text = format_certificate(results, "ATGGCTTAA", "ecoli")
        assert f"v{VERSION}" in text


# ---------------------------------------------------------------------------
# generate_certificate (full dataclass builder)
# ---------------------------------------------------------------------------

class TestGenerateCertificate:
    """``generate_certificate`` builds a serializable Certificate object."""

    def _to_type_results(self, predicate_results: list[PredicateResult]) -> list[TypeCheckResult]:
        """Convert PredicateResult -> TypeCheckResult for generate_certificate."""
        out = []
        for pr in predicate_results:
            v = pr.verdict if pr.verdict is not None else (
                Verdict.PASS if pr.passed else Verdict.FAIL
            )
            out.append(TypeCheckResult(predicate=pr.predicate, verdict=v))
        return out

    def _default_input_params(self, sequence: str) -> dict:
        return {
            "organism": "Homo_sapiens",
            "gc_lo": 0.30,
            "gc_hi": 0.70,
            "cai_threshold": 0.5,
            "enzymes": ["EcoRI", "BamHI", "XhoI", "HindIII", "NotI"],
            "exon_boundaries": [(0, len(sequence))],
        }

    def test_generate_certificate_returns_certificate_dataclass(self):
        seq = "ATGGTGCACCTGACCCCTGAGGAGAAG"
        type_results = self._to_type_results(_all_pass_results())
        cert = generate_certificate(
            seq, type_results, self._default_input_params(seq),
        )
        from biocompiler.shared.types import Certificate
        assert isinstance(cert, Certificate)

    def test_generate_certificate_has_required_fields(self):
        seq = "ATGGTGCACCTGACCCCTGAGGAGAAG"
        type_results = self._to_type_results(_all_pass_results())
        cert = generate_certificate(
            seq, type_results, self._default_input_params(seq),
        )
        # Certificate dataclass required fields.
        assert cert.version
        assert cert.design_id
        assert cert.sequence == seq
        assert isinstance(cert.types, list) and len(cert.types) == len(type_results)
        assert isinstance(cert.provenance, dict)
        # Provenance must record tool, version, timestamp, input_hash.
        for key in ("tool", "version", "timestamp", "input_hash"):
            assert key in cert.provenance, f"provenance missing key: {key}"

    def test_generate_certificate_embeds_verification_parameters(self):
        """The certificate must embed ALL parameters needed for
        independent re-verification (organism, gc_lo, gc_hi, cai_threshold,
        enzymes, exon_boundaries)."""
        seq = "ATGGTGCACCTGACCCCTGAGGAGAAG"
        type_results = self._to_type_results(_all_pass_results())
        cert = generate_certificate(
            seq, type_results, self._default_input_params(seq),
        )
        params = cert.provenance["parameters"]
        for key in ("organism", "gc_lo", "gc_hi", "cai_threshold", "enzymes",
                    "exon_boundaries", "cryptic_splice_threshold"):
            assert key in params, f"parameters missing key: {key}"

    def test_generate_certificate_records_solver_backend(self):
        """Provenance must record which solver backend was used so the
        certificate is reproducible."""
        seq = "ATGGTGCACCTGACCCCTGAGGAGAAG"
        type_results = self._to_type_results(_all_pass_results())
        cert = generate_certificate(
            seq, type_results, self._default_input_params(seq),
            solver_backend="greedy",
        )
        assert cert.provenance["solver_backend"] == "greedy"
        assert cert.provenance["solver_config"] == {}

    def test_generate_certificate_graduated_mode_does_not_raise(self):
        """GRADUATED mode (default) must NOT raise even if some
        predicates fail; the certificate documents the failures."""
        seq = "ATGGTGCACCTGACCCCTGAGGAGAAG"
        results = _all_pass_results()
        results[0] = PredicateResult(
            "NoStopCodons", False, verdict=Verdict.FAIL,
            details="in-frame stop found",
        )
        type_results = self._to_type_results(results)
        # Should NOT raise.
        cert = generate_certificate(
            seq, type_results, self._default_input_params(seq),
        )
        assert cert.provenance["overall_status"].startswith("PARTIAL_")

    def test_generate_certificate_strict_mode_raises_on_failure(self):
        """STRICT mode (require_all_pass=True) must raise
        CertificateGenerationError when any predicate fails."""
        from biocompiler.shared.exceptions import CertificateGenerationError

        seq = "ATGGTGCACCTGACCCCTGAGGAGAAG"
        results = _all_pass_results()
        results[0] = PredicateResult(
            "NoStopCodons", False, verdict=Verdict.FAIL,
            details="in-frame stop found",
        )
        type_results = self._to_type_results(results)
        with pytest.raises(CertificateGenerationError):
            generate_certificate(
                seq, type_results, self._default_input_params(seq),
                require_all_pass=True,
            )

    def test_generate_certificate_rejects_empty_sequence(self):
        type_results = self._to_type_results(_all_pass_results())
        with pytest.raises(ValueError, match="empty"):
            generate_certificate("", type_results, self._default_input_params(""))

    def test_generate_certificate_rejects_empty_type_results(self):
        seq = "ATGGTGCACCTGACCCCTGAGGAGAAG"
        with pytest.raises(ValueError, match="empty"):
            generate_certificate(seq, [], self._default_input_params(seq))

    def test_certificate_to_dict_roundtrip(self):
        """Certificate.to_dict() must produce a JSON-serializable dict
        that can be re-validated."""
        import json as _json
        seq = "ATGGTGCACCTGACCCCTGAGGAGAAG"
        type_results = self._to_type_results(_all_pass_results())
        cert = generate_certificate(
            seq, type_results, self._default_input_params(seq),
        )
        d = cert.to_dict()
        # Must be JSON-serializable (no datetime, no enums by value).
        _json.dumps(d)
        # Required top-level keys per _CERT_REQUIRED_KEYS.
        for key in ("version", "design_id", "sequence", "types", "provenance"):
            assert key in d


# ---------------------------------------------------------------------------
# verify_certificate (in-package independent re-verifier)
# ---------------------------------------------------------------------------

class TestVerifyCertificate:
    """``verify_certificate`` independently re-checks a certificate dict."""

    def _good_cert_dict(self) -> dict:
        seq = "ATGGTGCACCTGACCCCTGAGGAGAAG"
        type_results = [
            TypeCheckResult(predicate="NoStopCodons", verdict=Verdict.PASS),
            TypeCheckResult(predicate="ValidCodingSeq", verdict=Verdict.PASS),
        ]
        cert = generate_certificate(
            seq, type_results,
            {
                "organism": "Homo_sapiens", "gc_lo": 0.30, "gc_hi": 0.70,
                "cai_threshold": 0.5, "enzymes": [], "exon_boundaries": [(0, len(seq))],
            },
        )
        return cert.to_dict()

    def test_verify_certificate_returns_status_and_failures(self):
        d = self._good_cert_dict()
        status, failures = verify_certificate(d)
        assert status in ("VERIFIED", "REJECTED")
        assert isinstance(failures, list)

    def test_verify_certificate_rejects_missing_keys(self):
        """A cert dict missing required top-level keys is REJECTED."""
        status, failures = verify_certificate({})
        assert status == "REJECTED"
        assert any("missing" in f.lower() for f in failures)

    def test_verify_certificate_detects_hash_tamper(self):
        """If the design_id is mutated, verification must REJECT."""
        d = self._good_cert_dict()
        d["design_id"] = "0" * 64
        status, failures = verify_certificate(d)
        assert status == "REJECTED"
        assert any("design_id" in f.lower() and "mismatch" in f.lower()
                   for f in failures)


# ---------------------------------------------------------------------------
# End-to-end: optimizer -> certificate tier
# ---------------------------------------------------------------------------

class TestEndToEndCertificateFromOptimizer:
    """Optimize a gene, compute its certificate, check the tier."""

    @pytest.mark.parametrize("protein,organism", [
        ("MVHLTPEEK", "human"),
        ("MVSKGE", "human"),
        ("MFSF", "human"),
    ])
    def test_optimizer_produces_valid_certificate_tier(self, protein, organism):
        """Every optimizer run must produce a predicate_results list
        from which ``compute_certificate`` returns one of GOLD/SILVER/
        BRONZE (never None, never an unknown value)."""
        result = optimize_sequence(protein, organism=organism, strict_mode=False)
        assert result.predicate_results, "Optimizer produced no predicate results"
        tier = compute_certificate(result.predicate_results)
        assert tier in (CertLevel.GOLD, CertLevel.SILVER, CertLevel.BRONZE), (
            f"Unknown tier {tier!r} for protein {protein!r}"
        )

    def test_optimizer_certificate_text_matches_tier(self):
        """The optimizer's own ``certificate_text`` field must mention
        the same tier that ``compute_certificate`` returns."""
        result = optimize_sequence("MVHLTPEEK", organism="human", strict_mode=False)
        tier = compute_certificate(result.predicate_results)
        assert result.certificate_text, "Optimizer did not produce certificate_text"
        assert tier.value in result.certificate_text, (
            f"Tier {tier.value} not found in certificate_text"
        )

    def test_optimizer_certificate_text_lists_all_three_tiers(self):
        """The optimizer's certificate_text footer must explain GOLD/
        SILVER/BRONZE so consumers can interpret the level."""
        result = optimize_sequence("MVHLTPEEK", organism="human", strict_mode=False)
        text = result.certificate_text
        assert "GOLD" in text
        assert "SILVER" in text
        assert "BRONZE" in text

    def test_optimizer_predicate_results_have_required_fields(self):
        """Every PredicateResult from the optimizer must have a predicate
        name, a passed flag, and a details string."""
        result = optimize_sequence("MVHLTPEEK", organism="human", strict_mode=False)
        for pr in result.predicate_results:
            assert pr.predicate, "PredicateResult missing predicate name"
            assert isinstance(pr.passed, bool), (
                f"PredicateResult {pr.predicate} has non-bool passed"
            )
            assert isinstance(pr.details, str)

    def test_uncertain_count_caps_tier_correctly(self):
        """The optimizer's reported ``uncertain_predicate_count`` must
        drive the certificate tier per the capping rules:
            0 UNCERTAIN -> GOLD or SILVER (depending on mutagenesis)
            1 UNCERTAIN -> SILVER (or BRONZE if also failing)
            >=2 UNCERTAIN -> BRONZE
        """
        result = optimize_sequence("MVHLTPEEK", organism="human", strict_mode=False)
        tier = compute_certificate(result.predicate_results)
        n_uncertain = result.uncertain_predicate_count
        n_failed = len(result.failed_predicates)
        if n_failed > 0 or n_uncertain >= 2:
            assert tier == CertLevel.BRONZE, (
                f"{n_uncertain} UNCERTAIN + {n_failed} FAIL should be BRONZE, "
                f"got {tier.value}"
            )
        elif n_uncertain == 1:
            assert tier == CertLevel.SILVER, (
                f"1 UNCERTAIN should be SILVER, got {tier.value}"
            )


# ---------------------------------------------------------------------------
# End-to-end: optimizer -> standalone verifier (cross-check)
# ---------------------------------------------------------------------------

class TestEndToEndCrossCheckWithStandaloneVerifier:
    """Cross-check that the in-package certificate tier is consistent
    with the standalone verifier's verdict on the same optimized
    sequence."""

    def _try_import_standalone(self):
        try:
            import standalone_verifier as sv  # type: ignore[import-not-found]
            return sv
        except ImportError:
            pytest.skip("standalone_verifier not importable from this path")

    def test_standalone_verifier_accepts_optimizer_sequence(self):
        """The standalone verifier (a separate, stdlib-only
        implementation) must independently confirm that the optimizer's
        14 core DNA-level predicates all PASS on the optimized sequence."""
        sv = self._try_import_standalone()
        result = optimize_sequence("MVHLTPEEK", organism="human", strict_mode=False)
        cert = {
            "design_id": "",
            "organism": "Homo_sapiens",
            "protein_sequence": result.protein,
            "dna_sequence": result.sequence,
            "original_protein": result.protein,
            "original_dna_sequence": result.sequence,
            "gc_range": [0.30, 0.70],
            "certificate_level": compute_certificate(
                result.predicate_results
            ).value,
            "predicates": {},
            "optimizer_version": "1.0.0",
        }
        cert["design_id"] = sv.compute_design_id(cert)
        results = sv.verify_certificate(cert)
        # All 14 core predicates must PASS.
        core = {
            n: r for n, r in results.items()
            if not n.startswith("_") and n not in sv.SLOT_PREDICATES
        }
        failed = {n: r for n, r in core.items() if r["verdict"] != "PASS"}
        assert not failed, (
            f"Standalone verifier flagged core predicates as not-PASS: {failed}"
        )
        # Hash integrity must also PASS.
        assert results["_hash_integrity"]["verdict"] == "PASS"


if __name__ == "__main__":
    import pytest
    sys.exit(pytest.main([__file__, "-v"]))
