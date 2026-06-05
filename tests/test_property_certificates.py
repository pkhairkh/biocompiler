"""
Property-based tests for BioCompiler Certificate Computation Consistency.

Verifies that the Python certificate implementation matches the Lean4 theorems
from proof/BioCompiler/Certificates.lean (previously lean/BioCompiler/Certificate.lean), specifically:

  1. gold_implies_all_optimization   — GOLD ⇒ every predicate passed via optimization alone
  2. silver_implies_all_passed       — SILVER ⇒ no predicate is unsatisfied
  3. bronze_implies_unsatisfied      — BRONZE ⇒ at least one predicate is unsatisfied
  4. Determinism                     — same input always yields same output
  5. Certificate round-trip          — to_dict → from_dict → to_dict is identity
  6. Certificate validation          — missing keys are caught
  7. verify_certificate              — re-evaluates predicates independently
  8. All-pass ⇒ GOLD, any-fail ⇒ BRONZE
  9. Monotonicity                    — adding failures can only lower the level
"""

import hashlib
from datetime import datetime, timezone
from hypothesis import given, assume, settings, HealthCheck
from hypothesis import strategies as st
from hypothesis.strategies import SearchStrategy
import pytest

from biocompiler.certificate import (
    compute_certificate,
    generate_certificate,
    verify_certificate,
    _validate_cert_structure,
    _CERT_REQUIRED_KEYS,
    _PROVENANCE_REQUIRED_KEYS,
)
from biocompiler.type_system import CertLevel, PredicateResult, SpliceVerdict
from biocompiler.types import Verdict, TypeCheckResult, Certificate


# ────────────────────────────────────────────────────────────
# Hypothesis Strategies
# ────────────────────────────────────────────────────────────

PREDICATE_NAMES = [
    "NoStopCodons", "NoCrypticSplice", "NoCpGIsland",
    "NoRestrictionSite", "NoGTDinucleotide", "ValidCodingSeq",
    "ConservationScore", "CodonOptimality",
]

# Satisfaction method determines both `passed` and `details`
# to match the Lean4 SatisfactionMethod inductive type.
SATISFACTION_METHODS = ("optimization", "mutagenesis", "unsatisfied")


def _details_for_method(method: str, predicate: str) -> str:
    """Generate a plausible details string for a given satisfaction method."""
    if method == "optimization":
        return f"No issues found for {predicate}"
    elif method == "mutagenesis":
        return f"Resolved via mutagenesis for {predicate}"
    elif method == "unavoidable":
        return f"Unavoidable constraint in {predicate}"
    else:  # unsatisfied
        return f"Unsatisfied: {predicate} check failed"


@st.composite
def satisfaction_method_strategy(draw):
    """Draw a satisfaction method that maps to Lean4 SatisfactionMethod."""
    return draw(st.sampled_from(SATISFACTION_METHODS))


@st.composite
def predicate_result_strategy(draw):
    """Generate a PredicateResult consistent with a satisfaction method.

    This mirrors the Lean4 PredicateRecord where method determines passed:
      - optimization  → passed=True, no mutagenesis/unavoidable in details
      - mutagenesis   → passed=True, "mutagenesis" in details
      - unavoidable   → passed=True, "unavoidable" in details  (also → SILVER)
      - unsatisfied   → passed=False
    """
    method = draw(st.sampled_from(SATISFACTION_METHODS + ("unavoidable",)))
    predicate = draw(st.sampled_from(PREDICATE_NAMES))
    verdict = draw(st.sampled_from(list(Verdict)))

    if method == "optimization":
        passed = True
        details = f"No issues found for {predicate}"
    elif method == "mutagenesis":
        passed = True
        details = f"Resolved via mutagenesis for {predicate}"
    elif method == "unavoidable":
        passed = True
        details = f"Unavoidable constraint in {predicate}"
    else:  # unsatisfied
        passed = False
        details = f"Unsatisfied: {predicate} check failed"

    # Only attach SpliceVerdict to NoCrypticSplice
    splice_verdict = None
    if predicate == "NoCrypticSplice":
        splice_verdict = draw(st.sampled_from(list(SpliceVerdict)))

    return PredicateResult(
        predicate=predicate,
        passed=passed,
        verdict=verdict,
        details=details,
        positions=[],
    )


@st.composite
def predicate_results_list_strategy(draw, min_size=1, max_size=20):
    """Generate a list of PredicateResults."""
    return draw(
        st.lists(
            predicate_result_strategy(),
            min_size=min_size,
            max_size=max_size,
        )
    )


# Strategy for DNA sequences
dna_bases = st.sampled_from("ACGT")
dna_sequence = st.text(dna_bases, min_size=9, max_size=300)


@st.composite
def type_check_result_strategy(draw):
    """Generate a TypeCheckResult for generate_certificate testing."""
    predicate = draw(st.sampled_from(PREDICATE_NAMES))
    verdict = draw(st.sampled_from(list(Verdict)))
    return TypeCheckResult(
        predicate=predicate,
        verdict=verdict,
        derivation=None,
    )


@st.composite
def certificate_dict_strategy(draw):
    """Generate a well-formed certificate dict for round-trip testing."""
    seq = draw(dna_sequence)
    seq_hash = hashlib.sha256(seq.encode()).hexdigest()
    version = draw(st.from_regex(r"[0-9]+\.[0-9]+\.[0-9]+", fullmatch=True))
    design_id = seq_hash

    n_types = draw(st.integers(min_value=1, max_value=8))
    types_list = []
    for _ in range(n_types):
        predicate = draw(st.sampled_from(PREDICATE_NAMES))
        verdict_val = draw(st.sampled_from([v.value for v in Verdict]))
        types_list.append({
            "predicate": predicate,
            "verdict": verdict_val,
            "derivation": None,
        })

    provenance = {
        "tool": "BioCompiler",
        "version": version,
        "timestamp": datetime.now(timezone.utc).isoformat(),
        "input_hash": seq_hash,
    }

    return {
        "version": version,
        "design_id": design_id,
        "sequence": seq,
        "types": types_list,
        "provenance": provenance,
    }


# ────────────────────────────────────────────────────────────
# Theorem 1: gold_implies_all_optimization
# ────────────────────────────────────────────────────────────

class TestGoldImpliesAllOptimization:
    """Corresponds to Lean4: gold_implies_all_optimization

    theorem gold_implies_all_optimization
        (results : List PredicateRecord)
        (hcert : computeCertificate results = CertLevel.GOLD)
        (r : PredicateRecord) (hr : r ∈ results) :
        r.method = .optimization

    In Python: if compute_certificate(results) == GOLD, then every
    result passed=True AND no "mutagenesis"/"unavoidable" in details.
    """

    @given(results=predicate_results_list_strategy())
    @settings(max_examples=200, suppress_health_check=[HealthCheck.too_slow])
    def test_gold_implies_all_optimization(self, results):
        """If certificate is GOLD, every predicate was satisfied by optimization alone."""
        cert = compute_certificate(results)
        if cert == CertLevel.GOLD:
            for r in results:
                assert r.passed is True, (
                    f"GOLD certificate but predicate {r.predicate} not passed"
                )
                assert "mutagenesis" not in r.details.lower(), (
                    f"GOLD certificate but predicate {r.predicate} used mutagenesis: {r.details}"
                )
                assert "unavoidable" not in r.details.lower(), (
                    f"GOLD certificate but predicate {r.predicate} has unavoidable constraint: {r.details}"
                )

    @given(results=predicate_results_list_strategy())
    @settings(max_examples=200, suppress_health_check=[HealthCheck.too_slow])
    def test_gold_means_all_passed_no_mutagenesis(self, results):
        """GOLD iff all passed and no mutagenesis/unavoidable details."""
        cert = compute_certificate(results)
        all_passed = all(r.passed for r in results)
        no_mutagenesis = all("mutagenesis" not in r.details.lower() for r in results)
        no_unavoidable = all("unavoidable" not in r.details.lower() for r in results)

        if all_passed and no_mutagenesis and no_unavoidable:
            assert cert == CertLevel.GOLD, (
                f"All passed with no mutagenesis/unavoidable but cert={cert.value}"
            )


# ────────────────────────────────────────────────────────────
# Theorem 2: silver_implies_all_passed
# ────────────────────────────────────────────────────────────

class TestSilverImpliesAllPassed:
    """Corresponds to Lean4: silver_implies_all_passed

    theorem silver_implies_all_passed
        (results : List PredicateRecord)
        (hcert : computeCertificate results = CertLevel.SILVER)
        (hwf : ∀ r ∈ results, r.method ≠ .unsatisfied → r.passed = true)
        (r : PredicateRecord) (hr : r ∈ results) :
        r.passed = true ∧ r.method ≠ .unsatisfied

    In Python: if SILVER, then every result passed=True.
    The well-formedness condition (hwf) in Lean maps to the Python invariant
    that any non-unsatisfied predicate should have passed=True, which is
    guaranteed by our strategy construction.
    """

    @given(results=predicate_results_list_strategy())
    @settings(max_examples=200, suppress_health_check=[HealthCheck.too_slow])
    def test_silver_implies_all_passed(self, results):
        """If certificate is SILVER, all predicates passed."""
        cert = compute_certificate(results)
        if cert == CertLevel.SILVER:
            for r in results:
                assert r.passed is True, (
                    f"SILVER certificate but predicate {r.predicate} not passed"
                )

    @given(results=predicate_results_list_strategy())
    @settings(max_examples=200, suppress_health_check=[HealthCheck.too_slow])
    def test_silver_implies_no_unsatisfied(self, results):
        """If certificate is SILVER, no predicate is unsatisfied (passed=False)."""
        cert = compute_certificate(results)
        if cert == CertLevel.SILVER:
            has_unsatisfied = any(not r.passed for r in results)
            assert not has_unsatisfied, (
                "SILVER certificate but at least one predicate is unsatisfied"
            )

    @given(results=predicate_results_list_strategy())
    @settings(max_examples=200, suppress_health_check=[HealthCheck.too_slow])
    def test_silver_has_mutagenesis_or_unavoidable(self, results):
        """SILVER iff all passed AND at least one has mutagenesis/unavoidable."""
        cert = compute_certificate(results)
        all_passed = all(r.passed for r in results)
        has_mutagenesis = any("mutagenesis" in r.details.lower() for r in results)
        has_unavoidable = any("unavoidable" in r.details.lower() for r in results)

        if all_passed and (has_mutagenesis or has_unavoidable):
            # Note: could also be BRONZE if some are unsatisfied,
            # but all_passed excludes that. So it must be SILVER.
            assert cert == CertLevel.SILVER, (
                f"All passed with mutagenesis/unavoidable but cert={cert.value}"
            )


# ────────────────────────────────────────────────────────────
# Theorem 3: bronze_implies_unsatisfied
# ────────────────────────────────────────────────────────────

class TestBronzeImpliesUnsatisfied:
    """Corresponds to Lean4: bronze_implies_unsatisfied

    theorem bronze_implies_unsatisfied
        (results : List PredicateRecord)
        (hcert : computeCertificate results = CertLevel.BRONZE) :
        ∃ r ∈ results, r.method = .unsatisfied

    In Python: if BRONZE, then at least one result has passed=False.
    """

    @given(results=predicate_results_list_strategy())
    @settings(max_examples=200, suppress_health_check=[HealthCheck.too_slow])
    def test_bronze_implies_at_least_one_unsatisfied(self, results):
        """If certificate is BRONZE, at least one predicate is unsatisfied."""
        cert = compute_certificate(results)
        if cert == CertLevel.BRONZE:
            unsatisfied = [r for r in results if not r.passed]
            assert len(unsatisfied) >= 1, (
                "BRONZE certificate but no unsatisfied predicates found"
            )

    @given(results=predicate_results_list_strategy())
    @settings(max_examples=200, suppress_health_check=[HealthCheck.too_slow])
    def test_any_unsatisfied_implies_bronze(self, results):
        """If any predicate is unsatisfied (passed=False), cert must be BRONZE."""
        has_unsatisfied = any(not r.passed for r in results)
        cert = compute_certificate(results)
        if has_unsatisfied:
            assert cert == CertLevel.BRONZE, (
                f"Has unsatisfied predicates but cert={cert.value}"
            )


# ────────────────────────────────────────────────────────────
# Theorem: Determinism
# ────────────────────────────────────────────────────────────

class TestDeterminism:
    """compute_certificate is a pure function: same input → same output."""

    @given(results=predicate_results_list_strategy())
    @settings(max_examples=200, suppress_health_check=[HealthCheck.too_slow])
    def test_compute_certificate_is_deterministic(self, results):
        """Calling compute_certificate twice with same input yields same result."""
        cert1 = compute_certificate(results)
        cert2 = compute_certificate(results)
        assert cert1 == cert2

    @given(results=predicate_results_list_strategy())
    @settings(max_examples=100, suppress_health_check=[HealthCheck.too_slow])
    def test_compute_certificate_deterministic_many_calls(self, results):
        """Calling compute_certificate 10 times always yields the same result."""
        certs = [compute_certificate(results) for _ in range(10)]
        assert all(c == certs[0] for c in certs)


# ────────────────────────────────────────────────────────────
# Theorem: Certificate Round-Trip
# ────────────────────────────────────────────────────────────

class TestCertificateRoundTrip:
    """Certificate.to_dict → Certificate.from_dict → to_dict is identity."""

    @given(cert_dict=certificate_dict_strategy())
    @settings(max_examples=100, suppress_health_check=[HealthCheck.too_slow])
    def test_round_trip_to_dict_from_dict(self, cert_dict):
        """to_dict → from_dict → to_dict yields identical dict."""
        cert = Certificate.from_dict(cert_dict)
        round_tripped = cert.to_dict()
        assert round_tripped == cert_dict

    @given(cert_dict=certificate_dict_strategy())
    @settings(max_examples=100, suppress_health_check=[HealthCheck.too_slow])
    def test_from_dict_preserves_all_fields(self, cert_dict):
        """from_dict preserves all required fields."""
        cert = Certificate.from_dict(cert_dict)
        assert cert.version == cert_dict["version"]
        assert cert.design_id == cert_dict["design_id"]
        assert cert.sequence == cert_dict["sequence"]
        assert cert.types == cert_dict["types"]
        assert cert.provenance == cert_dict["provenance"]


# ────────────────────────────────────────────────────────────
# Theorem: Certificate Validation Catches Missing Keys
# ────────────────────────────────────────────────────────────

class TestCertificateValidation:
    """Missing required keys are detected by _validate_cert_structure."""

    @given(
        missing_key=st.sampled_from(sorted(_CERT_REQUIRED_KEYS)),
    )
    @settings(max_examples=50, suppress_health_check=[HealthCheck.too_slow])
    def test_missing_cert_key_detected(self, missing_key):
        """Removing a required certificate key triggers validation failure."""
        cert_dict = {
            "version": "7.0.0",
            "design_id": "abc123",
            "sequence": "ATGGCTTAA",
            "types": [],
            "provenance": {
                "tool": "BioCompiler",
                "version": "7.0.0",
                "timestamp": "2026-01-01T00:00:00Z",
                "input_hash": "abc123",
            },
        }
        del cert_dict[missing_key]
        issues = _validate_cert_structure(cert_dict)
        assert len(issues) > 0
        assert missing_key in issues[0]

    @given(
        missing_prov_key=st.sampled_from(sorted(_PROVENANCE_REQUIRED_KEYS)),
    )
    @settings(max_examples=50, suppress_health_check=[HealthCheck.too_slow])
    def test_missing_provenance_key_detected(self, missing_prov_key):
        """Removing a required provenance key triggers validation failure."""
        cert_dict = {
            "version": "7.0.0",
            "design_id": "abc123",
            "sequence": "ATGGCTTAA",
            "types": [],
            "provenance": {
                "tool": "BioCompiler",
                "version": "7.0.0",
                "timestamp": "2026-01-01T00:00:00Z",
                "input_hash": "abc123",
            },
        }
        del cert_dict["provenance"][missing_prov_key]
        issues = _validate_cert_structure(cert_dict)
        assert len(issues) > 0

    def test_missing_predicate_key_in_type_entry(self):
        """Type entry missing 'predicate' or 'verdict' is caught."""
        cert_dict = {
            "version": "7.0.0",
            "design_id": "abc123",
            "sequence": "ATGGCTTAA",
            "types": [{"predicate": "NoStopCodons"}],  # missing verdict
            "provenance": {
                "tool": "BioCompiler",
                "version": "7.0.0",
                "timestamp": "2026-01-01T00:00:00Z",
                "input_hash": "abc123",
            },
        }
        issues = _validate_cert_structure(cert_dict)
        assert len(issues) > 0
        assert "verdict" in issues[0]

    def test_valid_certificate_has_no_issues(self):
        """A well-formed certificate passes validation."""
        cert_dict = {
            "version": "7.0.0",
            "design_id": "abc123",
            "sequence": "ATGGCTTAA",
            "types": [
                {"predicate": "NoStopCodons", "verdict": "PASS"},
            ],
            "provenance": {
                "tool": "BioCompiler",
                "version": "7.0.0",
                "timestamp": "2026-01-01T00:00:00Z",
                "input_hash": "abc123",
            },
        }
        issues = _validate_cert_structure(cert_dict)
        assert issues == []

    def test_from_dict_raises_on_missing_keys(self):
        """Certificate.from_dict raises ValueError when required keys are missing."""
        with pytest.raises(ValueError, match="missing keys"):
            Certificate.from_dict({"version": "7.0.0"})  # Missing most keys


# ────────────────────────────────────────────────────────────
# Theorem: All-pass ⇒ GOLD
# ────────────────────────────────────────────────────────────

class TestAllPassImpliesGold:
    """If all predicates pass (with no mutagenesis/unavoidable), cert is GOLD."""

    @given(
        n=st.integers(min_value=1, max_value=12),
    )
    @settings(max_examples=100, suppress_health_check=[HealthCheck.too_slow])
    def test_all_pass_yields_gold(self, n):
        """All predicates passed with no mutagenesis/unavoidable → GOLD."""
        results = [
            PredicateResult(
                predicate=f"Predicate_{i}",
                passed=True,
                verdict=Verdict.PASS,
                details="All good",
                positions=[],
            )
            for i in range(n)
        ]
        cert = compute_certificate(results)
        assert cert == CertLevel.GOLD

    @given(results=predicate_results_list_strategy())
    @settings(max_examples=200, suppress_health_check=[HealthCheck.too_slow])
    def test_gold_necessary_and_sufficient(self, results):
        """GOLD ⇔ all passed AND no mutagenesis AND no unavoidable."""
        cert = compute_certificate(results)
        all_passed = all(r.passed for r in results)
        no_mutagenesis = all("mutagenesis" not in r.details.lower() for r in results)
        no_unavoidable = all("unavoidable" not in r.details.lower() for r in results)

        is_gold_equivalent = all_passed and no_mutagenesis and no_unavoidable
        assert (cert == CertLevel.GOLD) == is_gold_equivalent, (
            f"GOLD equivalence mismatch: cert={cert.value}, "
            f"all_passed={all_passed}, no_mut={no_mutagenesis}, no_unavoid={no_unavoidable}"
        )


# ────────────────────────────────────────────────────────────
# Theorem: Any-fail ⇒ BRONZE
# ────────────────────────────────────────────────────────────

class TestAnyFailImpliesBronze:
    """If any predicate fails (passed=False), cert must be BRONZE."""

    @given(
        n_pass=st.integers(min_value=0, max_value=10),
        n_fail=st.integers(min_value=1, max_value=5),
    )
    @settings(max_examples=100, suppress_health_check=[HealthCheck.too_slow])
    def test_any_fail_yields_bronze(self, n_pass, n_fail):
        """Any predicate with passed=False → BRONZE."""
        results = []
        for i in range(n_pass):
            results.append(PredicateResult(
                predicate=f"PassPred_{i}",
                passed=True,
                verdict=Verdict.PASS,
                details="All good",
                positions=[],
            ))
        for i in range(n_fail):
            results.append(PredicateResult(
                predicate=f"FailPred_{i}",
                passed=False,
                verdict=Verdict.FAIL,
                details="Unsatisfied: check failed",
                positions=[42],
            ))
        cert = compute_certificate(results)
        assert cert == CertLevel.BRONZE


# ────────────────────────────────────────────────────────────
# Theorem: Monotonicity of Certificate Level
# ────────────────────────────────────────────────────────────

class TestMonotonicity:
    """Certificate level is monotonic in predicate satisfaction.

    Adding a failure can only lower the level. Removing a failure can
    only raise the level.  Order: GOLD > SILVER > BRONZE.
    """

    def _level_rank(self, level: CertLevel) -> int:
        """Numeric ranking: GOLD=2, SILVER=1, BRONZE=0."""
        return {CertLevel.GOLD: 2, CertLevel.SILVER: 1, CertLevel.BRONZE: 0}[level]

    @given(results=predicate_results_list_strategy())
    @settings(max_examples=200, suppress_health_check=[HealthCheck.too_slow])
    def test_adding_failure_lowers_or_maintains_level(self, results):
        """Adding a failing predicate cannot raise the certificate level."""
        original_cert = compute_certificate(results)

        # Add a failing predicate
        new_results = results + [
            PredicateResult(
                predicate="ExtraFail",
                passed=False,
                verdict=Verdict.FAIL,
                details="Unsatisfied: added failure",
                positions=[],
            )
        ]
        new_cert = compute_certificate(new_results)

        assert self._level_rank(new_cert) <= self._level_rank(original_cert), (
            f"Adding a failure raised level: {original_cert.value} → {new_cert.value}"
        )

    @given(results=predicate_results_list_strategy())
    @settings(max_examples=200, suppress_health_check=[HealthCheck.too_slow])
    def test_adding_passing_predicate_cannot_lower_level(self, results):
        """Adding a pure-optimization passing predicate cannot lower the level."""
        original_cert = compute_certificate(results)

        # Add a pure optimization-pass predicate
        new_results = results + [
            PredicateResult(
                predicate="ExtraPass",
                passed=True,
                verdict=Verdict.PASS,
                details="No issues found",
                positions=[],
            )
        ]
        new_cert = compute_certificate(new_results)

        assert self._level_rank(new_cert) >= self._level_rank(original_cert), (
            f"Adding a pass lowered level: {original_cert.value} → {new_cert.value}"
        )

    @given(results=predicate_results_list_strategy())
    @settings(max_examples=200, suppress_health_check=[HealthCheck.too_slow])
    def test_removing_failures_can_only_raise_level(self, results):
        """Removing all failing predicates can only raise or maintain the level."""
        original_cert = compute_certificate(results)

        # Remove all failing predicates
        passing_only = [r for r in results if r.passed]
        if not passing_only:
            # Can't have empty results, skip
            assume(len(passing_only) > 0)

        new_cert = compute_certificate(passing_only)
        assert self._level_rank(new_cert) >= self._level_rank(original_cert), (
            f"Removing failures lowered level: {original_cert.value} → {new_cert.value}"
        )


# ────────────────────────────────────────────────────────────
# Theorem: Exhaustiveness — Every result list yields exactly one level
# ────────────────────────────────────────────────────────────

class TestExhaustiveness:
    """Every non-empty list of PredicateResults maps to exactly one CertLevel."""

    @given(results=predicate_results_list_strategy())
    @settings(max_examples=200, suppress_health_check=[HealthCheck.too_slow])
    def test_compute_certificate_returns_valid_level(self, results):
        """compute_certificate always returns GOLD, SILVER, or BRONZE."""
        cert = compute_certificate(results)
        assert cert in (CertLevel.GOLD, CertLevel.SILVER, CertLevel.BRONZE)

    @given(results=predicate_results_list_strategy())
    @settings(max_examples=200, suppress_health_check=[HealthCheck.too_slow])
    def test_certificate_levels_are_mutually_exclusive(self, results):
        """The certificate level is exactly one of GOLD/SILVER/BRONZE."""
        cert = compute_certificate(results)
        levels = {CertLevel.GOLD, CertLevel.SILVER, CertLevel.BRONZE}
        assert cert in levels
        # It's exactly one (enum guarantees this, but we verify)
        matching = sum(1 for lv in levels if cert == lv)
        assert matching == 1


# ────────────────────────────────────────────────────────────
# Theorem: generate_certificate produces valid structure
# ────────────────────────────────────────────────────────────

class TestGenerateCertificate:
    """generate_certificate produces structurally valid certificates."""

    @given(
        seq=dna_sequence,
        type_results=st.lists(type_check_result_strategy(), min_size=1, max_size=8),
    )
    @settings(max_examples=50, suppress_health_check=[HealthCheck.too_slow])
    def test_generated_certificate_validates(self, seq, type_results):
        """A generated certificate passes structural validation."""
        cert = generate_certificate(
            sequence=seq,
            type_results=type_results,
            input_params={},
        )
        cert_dict = cert.to_dict()
        issues = _validate_cert_structure(cert_dict)
        assert issues == [], f"Generated certificate has structural issues: {issues}"

    @given(
        seq=dna_sequence,
        type_results=st.lists(type_check_result_strategy(), min_size=1, max_size=8),
    )
    @settings(max_examples=50, suppress_health_check=[HealthCheck.too_slow])
    def test_generated_certificate_hash_matches(self, seq, type_results):
        """The design_id in a generated certificate matches SHA-256 of sequence."""
        cert = generate_certificate(
            sequence=seq,
            type_results=type_results,
            input_params={},
        )
        expected_hash = hashlib.sha256(seq.encode()).hexdigest()
        assert cert.design_id == expected_hash

    @given(
        seq=dna_sequence,
        type_results=st.lists(type_check_result_strategy(), min_size=1, max_size=8),
    )
    @settings(max_examples=50, suppress_health_check=[HealthCheck.too_slow])
    def test_generated_certificate_round_trips(self, seq, type_results):
        """Generated certificate survives to_dict → from_dict round-trip."""
        cert = generate_certificate(
            sequence=seq,
            type_results=type_results,
            input_params={},
        )
        cert_dict = cert.to_dict()
        restored = Certificate.from_dict(cert_dict)
        assert restored.to_dict() == cert_dict


# ────────────────────────────────────────────────────────────
# Theorem: verify_certificate catches hash tampering
# ────────────────────────────────────────────────────────────

class TestVerifyCertificate:
    """verify_certificate detects structural and integrity issues."""

    def test_tampered_sequence_rejected(self):
        """Changing the sequence after signing should cause verification failure."""
        seq = "ATGGCTTAA"
        cert = generate_certificate(
            sequence=seq,
            type_results=[TypeCheckResult("NoStopCodons", Verdict.PASS)],
            input_params={},
        )
        cert_dict = cert.to_dict()
        # Tamper with the sequence but keep the old design_id
        cert_dict["sequence"] = "ATGGCTCCC"
        status, reasons = verify_certificate(cert_dict)
        assert status == "REJECTED"

    def test_missing_cert_key_rejected(self):
        """A certificate missing required keys is rejected."""
        status, reasons = verify_certificate({"version": "7.0.0"})
        assert status == "REJECTED"
        assert len(reasons) > 0

    def test_empty_types_not_rejected_structurally(self):
        """An empty types list passes structural validation (just no predicates checked)."""
        cert_dict = {
            "version": "7.0.0",
            "design_id": hashlib.sha256(b"ATGGCTTAA").hexdigest(),
            "sequence": "ATGGCTTAA",
            "types": [],
            "provenance": {
                "tool": "BioCompiler",
                "version": "7.0.0",
                "timestamp": "2026-01-01T00:00:00Z",
                "input_hash": hashlib.sha256(b"ATGGCTTAA").hexdigest(),
            },
        }
        # Structural check should pass (empty types is valid)
        issues = _validate_cert_structure(cert_dict)
        assert issues == []


# ────────────────────────────────────────────────────────────
# Theorem: Edge cases
# ────────────────────────────────────────────────────────────

class TestEdgeCases:
    """Edge cases and boundary conditions."""

    def test_single_passing_predicate(self):
        """Single passing predicate → GOLD."""
        results = [PredicateResult("NoStopCodons", True, details="OK")]
        assert compute_certificate(results) == CertLevel.GOLD

    def test_single_failing_predicate(self):
        """Single failing predicate → BRONZE."""
        results = [PredicateResult("NoStopCodons", False, details="Unsatisfied")]
        assert compute_certificate(results) == CertLevel.BRONZE

    def test_single_mutagenesis_predicate(self):
        """Single mutagenesis-resolved predicate → SILVER."""
        results = [PredicateResult("NoGTDinucleotide", True, details="mutagenesis applied")]
        assert compute_certificate(results) == CertLevel.SILVER

    def test_single_unavoidable_predicate(self):
        """Single unavoidable-resolved predicate → SILVER."""
        results = [PredicateResult("NoGTDinucleotide", True, details="unavoidable GT dinucleotide")]
        assert compute_certificate(results) == CertLevel.SILVER

    @given(
        n=st.integers(min_value=1, max_value=50),
    )
    @settings(max_examples=50, suppress_health_check=[HealthCheck.too_slow])
    def test_all_passing_large_list(self, n):
        """Large lists of all-passing predicates → GOLD."""
        results = [
            PredicateResult(f"P{i}", True, details="OK") for i in range(n)
        ]
        assert compute_certificate(results) == CertLevel.GOLD

    @given(
        n=st.integers(min_value=1, max_value=50),
    )
    @settings(max_examples=50, suppress_health_check=[HealthCheck.too_slow])
    def test_all_failing_large_list(self, n):
        """Large lists of all-failing predicates → BRONZE."""
        results = [
            PredicateResult(f"P{i}", False, details="Unsatisfied") for i in range(n)
        ]
        assert compute_certificate(results) == CertLevel.BRONZE
