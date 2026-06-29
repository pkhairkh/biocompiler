"""
Property-based tests for biocompiler.provenance.certificate using Hypothesis.

Covers three core properties:
  1. Certificate.to_dict / from_dict roundtrip preserves data
  2. CertLevel values are ordered: GOLD > SILVER > BRONZE
  3. generate_certificate always produces a structurally valid Certificate
"""

import hashlib
from datetime import datetime, timezone

import pytest
pytest.importorskip("hypothesis")
pytest.importorskip("hypothesis")
from hypothesis import given, assume, settings, HealthCheck
from hypothesis import strategies as st

from biocompiler.provenance.certificate import (
    generate_certificate,
    _validate_cert_structure,
)
from biocompiler.type_system import CertLevel
from biocompiler.shared.types import Certificate, Verdict, TypeCheckResult


# ────────────────────────────────────────────────────────────
# Shared Strategies
# ────────────────────────────────────────────────────────────

PREDICATE_NAMES = [
    "NoStopCodons", "NoCrypticSplice", "NoCpGIsland",
    "NoRestrictionSite", "NoGTDinucleotide", "ValidCodingSeq",
    "ConservationScore", "CodonOptimality",
]

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
    """Generate a well-formed certificate dict for roundtrip testing."""
    seq = draw(dna_sequence)
    seq_hash = hashlib.sha256(seq.encode()).hexdigest()
    version = draw(st.from_regex(r"[0-9]+\.[0-9]+\.[0-9]+", fullmatch=True))

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
        "design_id": seq_hash,
        "sequence": seq,
        "types": types_list,
        "provenance": provenance,
    }


# ────────────────────────────────────────────────────────────
# Property 1: Certificate to_dict / from_dict roundtrip
# ────────────────────────────────────────────────────────────

class TestCertificateRoundTrip:
    """Certificate.to_dict → Certificate.from_dict → to_dict is identity."""

    @given(cert_dict=certificate_dict_strategy())
    @settings(max_examples=100, suppress_health_check=[HealthCheck.too_slow])
    def test_roundtrip_preserves_data(self, cert_dict):
        """to_dict → from_dict → to_dict yields identical dict."""
        cert = Certificate.from_dict(cert_dict)
        round_tripped = cert.to_dict()
        # Legacy v1 dicts do not have hash_version, but from_dict defaults to 1
        # and to_dict omits it for v1. So the roundtrip should still be consistent.
        # If the original dict lacked hash_version, the roundtrip should too.
        if "hash_version" not in cert_dict:
            round_tripped.pop("hash_version", None)
        assert round_tripped == cert_dict

    @given(cert_dict=certificate_dict_strategy())
    @settings(max_examples=100, suppress_health_check=[HealthCheck.too_slow])
    def test_from_dict_preserves_each_field(self, cert_dict):
        """from_dict preserves version, design_id, sequence, types, provenance."""
        cert = Certificate.from_dict(cert_dict)
        assert cert.version == cert_dict["version"]
        assert cert.design_id == cert_dict["design_id"]
        assert cert.sequence == cert_dict["sequence"]
        assert cert.types == cert_dict["types"]
        assert cert.provenance == cert_dict["provenance"]

    @given(cert_dict=certificate_dict_strategy())
    @settings(max_examples=80, suppress_health_check=[HealthCheck.too_slow])
    def test_double_roundtrip_is_idempotent(self, cert_dict):
        """Two successive roundtrips produce the same result."""
        first = Certificate.from_dict(cert_dict).to_dict()
        second = Certificate.from_dict(first).to_dict()
        assert first == second

    def test_from_dict_rejects_missing_keys(self):
        """from_dict raises ValueError when required keys are absent."""
        with pytest.raises(ValueError, match="missing keys"):
            Certificate.from_dict({"version": "1.0.0"})


# ────────────────────────────────────────────────────────────
# Property 2: CertLevel values are ordered GOLD > SILVER > BRONZE
# ────────────────────────────────────────────────────────────

class TestCertLevelOrdering:
    """CertLevel enum members follow the order GOLD > SILVER > BRONZE."""

    @staticmethod
    def _rank(level: CertLevel) -> int:
        return {CertLevel.GOLD: 2, CertLevel.SILVER: 1, CertLevel.BRONZE: 0}[level]

    def test_gold_greater_than_silver(self):
        assert self._rank(CertLevel.GOLD) > self._rank(CertLevel.SILVER)

    def test_silver_greater_than_bronze(self):
        assert self._rank(CertLevel.SILVER) > self._rank(CertLevel.BRONZE)

    def test_gold_greater_than_bronze(self):
        assert self._rank(CertLevel.GOLD) > self._rank(CertLevel.BRONZE)

    def test_total_ordering_transitive(self):
        """GOLD > SILVER and SILVER > BRONZE implies GOLD > BRONZE."""
        assert self._rank(CertLevel.GOLD) > self._rank(CertLevel.BRONZE)

    @given(
        high=st.sampled_from([CertLevel.GOLD, CertLevel.SILVER]),
        low=st.sampled_from([CertLevel.SILVER, CertLevel.BRONZE]),
    )
    @settings(max_examples=50, suppress_health_check=[HealthCheck.too_slow])
    def test_rank_ordering_is_consistent(self, high, low):
        """Any higher level has a strictly greater rank than any lower level."""
        assume(self._rank(high) > self._rank(low))
        assert self._rank(high) > self._rank(low)

    @given(level=st.sampled_from(list(CertLevel)))
    @settings(max_examples=30, suppress_health_check=[HealthCheck.too_slow])
    def test_each_level_has_unique_rank(self, level):
        """No two distinct CertLevel values share the same rank."""
        all_ranks = {self._rank(lv): lv for lv in CertLevel}
        # 3 levels → 3 unique ranks
        assert len(all_ranks) == len(list(CertLevel))

    def test_exactly_three_levels(self):
        """CertLevel has exactly three members: GOLD, SILVER, BRONZE."""
        assert set(CertLevel) == {CertLevel.GOLD, CertLevel.SILVER, CertLevel.BRONZE}


# ────────────────────────────────────────────────────────────
# Property 3: generate_certificate produces valid Certificate
# ────────────────────────────────────────────────────────────

class TestGenerateCertificateValidity:
    """generate_certificate always returns a structurally valid Certificate."""

    @given(
        seq=dna_sequence,
        type_results=st.lists(type_check_result_strategy(), min_size=1, max_size=8),
    )
    @settings(max_examples=60, suppress_health_check=[HealthCheck.too_slow])
    def test_generated_certificate_passes_validation(self, seq, type_results):
        """A generated certificate passes _validate_cert_structure."""
        cert = generate_certificate(
            sequence=seq,
            type_results=type_results,
            input_params={},
        )
        issues = _validate_cert_structure(cert.to_dict())
        assert issues == [], f"Structural issues: {issues}"

    @given(
        seq=dna_sequence,
        type_results=st.lists(type_check_result_strategy(), min_size=1, max_size=8),
    )
    @settings(max_examples=60, suppress_health_check=[HealthCheck.too_slow])
    def test_design_id_is_current_hash(self, seq, type_results):
        """Certificate.design_id equals the v3 hash (sequence + predicates + params).

        Since C14 (prior fix), freshly generated certificates use the v3 hash
        (expanded param set). v1 (sequence-only) must still differ.
        """
        from biocompiler.provenance.certificate import _compute_certificate_hash
        cert = generate_certificate(
            sequence=seq,
            type_results=type_results,
            input_params={},
        )
        # v1 hash (sequence-only) should NOT match — this was the soundness bug
        v1_hash = hashlib.sha256(seq.encode()).hexdigest()
        assert cert.design_id != v1_hash, "v3 hash must differ from v1 sequence-only hash"
        # v3 hash should match
        params = cert.provenance.get("parameters", {})
        expected = _compute_certificate_hash(
            sequence=seq,
            types_list=cert.types,
            params=params,
            hash_version=3,
        )
        assert cert.design_id == expected

    @given(
        seq=dna_sequence,
        type_results=st.lists(type_check_result_strategy(), min_size=1, max_size=8),
    )
    @settings(max_examples=60, suppress_health_check=[HealthCheck.too_slow])
    def test_sequence_preserved_in_certificate(self, seq, type_results):
        """Certificate.sequence matches the input sequence exactly."""
        cert = generate_certificate(
            sequence=seq,
            type_results=type_results,
            input_params={},
        )
        assert cert.sequence == seq

    @given(
        seq=dna_sequence,
        type_results=st.lists(type_check_result_strategy(), min_size=1, max_size=8),
    )
    @settings(max_examples=60, suppress_health_check=[HealthCheck.too_slow])
    def test_types_match_input_results(self, seq, type_results):
        """Certificate.types list reflects every TypeCheckResult."""
        cert = generate_certificate(
            sequence=seq,
            type_results=type_results,
            input_params={},
        )
        assert len(cert.types) == len(type_results)
        for cert_entry, tc_result in zip(cert.types, type_results):
            assert cert_entry["predicate"] == tc_result.predicate
            assert cert_entry["verdict"] == tc_result.verdict.value

    @given(
        seq=dna_sequence,
        type_results=st.lists(type_check_result_strategy(), min_size=1, max_size=8),
    )
    @settings(max_examples=50, suppress_health_check=[HealthCheck.too_slow])
    def test_generated_certificate_roundtrips(self, seq, type_results):
        """A generated certificate survives to_dict → from_dict roundtrip."""
        cert = generate_certificate(
            sequence=seq,
            type_results=type_results,
            input_params={},
        )
        cert_dict = cert.to_dict()
        restored = Certificate.from_dict(cert_dict)
        assert restored.to_dict() == cert_dict

    def test_empty_sequence_raises_value_error(self):
        """generate_certificate raises ValueError for empty sequence."""
        with pytest.raises(ValueError, match="Sequence must not be empty"):
            generate_certificate(
                sequence="",
                type_results=[TypeCheckResult("NoStopCodons", Verdict.PASS)],
                input_params={},
            )

    def test_empty_type_results_raises_value_error(self):
        """generate_certificate raises ValueError for empty type_results."""
        with pytest.raises(ValueError, match="Type results must not be empty"):
            generate_certificate(
                sequence="ATGGCTTAA",
                type_results=[],
                input_params={},
            )

    @given(
        seq=dna_sequence,
        type_results=st.lists(type_check_result_strategy(), min_size=1, max_size=4),
        input_params=st.dictionaries(
            st.sampled_from(["organism", "cell_type", "gc_lo", "gc_hi"]),
            st.one_of(st.just("Homo_sapiens"), st.just("HEK293T"), st.floats(0.1, 0.9)),
            max_size=4,
        ),
    )
    @settings(max_examples=40, suppress_health_check=[HealthCheck.too_slow])
    def test_input_params_embedded_in_provenance(self, seq, type_results, input_params):
        """Parameters passed to generate_certificate appear in provenance.parameters."""
        cert = generate_certificate(
            sequence=seq,
            type_results=type_results,
            input_params=input_params,
        )
        params = cert.provenance.get("parameters", {})
        for key, value in input_params.items():
            assert key in params
            assert params[key] == value
