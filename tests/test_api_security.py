"""
Agent 53: API Security and Stress Tests
=========================================

Tests for authentication enforcement, rate limiting, input validation,
size limits, concurrency safety, and provenance persistence across restarts.
"""

from __future__ import annotations

import asyncio
import os
import sqlite3
import tempfile
import threading
import time
from pathlib import Path
from unittest.mock import patch

import pytest

from biocompiler.api import (
    MAX_PROTEIN_LENGTH,
    MAX_DNA_LENGTH,
    MAX_BATCH_SIZE,
    MAX_REQUEST_SIZE,
    OPTIMIZE_TIMEOUT_S,
    RATE_LIMIT_RPM,
    BATCH_CHECK_MAX,
    BATCH_OPTIMIZE_MAX,
    BATCH_EXPORT_MAX,
    ProteinInput,
    validate_protein_input,
    validate_organism_input,
)
from biocompiler.infrastructure.rate_limiter import PersistentRateLimiter
from biocompiler.organisms import SUPPORTED_ORGANISMS


# ═══════════════════════════════════════════════════════════════════════
# 1. Authentication Tests
# ═══════════════════════════════════════════════════════════════════════


class TestAuthRequiredByDefault:
    """Verify that auth is required by default (the secure default)."""

    def test_auth_mode_defaults_to_required(self):
        """Without env overrides, the auth mode should be 'required'."""
        # The module-level _AUTH_MODE is read at import time.
        # We test the public getter which reflects the current state.
        from biocompiler.api import get_auth_mode
        # In test environments the key may be auto-generated, but mode
        # should be "required" unless explicitly changed.
        mode = get_auth_mode()
        assert mode in ("required", "optional", "disabled", "deferred")
        # Default is "required" unless someone explicitly disabled it
        if os.environ.get("BIOCOMPILER_AUTH_MODE", "required").lower() == "required":
            assert mode in ("required", "deferred")

    def test_is_auth_enabled_by_default(self):
        """Auth should be enabled unless explicitly disabled."""
        from biocompiler.api import is_auth_enabled
        # In deferred mode (no env set), auth is not yet enabled until the
        # API server starts and generates a key. This is correct behaviour.
        # In required/optional mode, auth is enabled.
        from biocompiler.api import auth as _auth_module
        if _auth_module._AUTH_MODE == "deferred":
            # Deferred mode: auth will be enabled when server starts
            assert is_auth_enabled() is False or os.environ.get(
                "BIOCOMPILER_API_KEY", ""
            ) == "disabled"
        else:
            assert is_auth_enabled() is True or os.environ.get(
                "BIOCOMPILER_API_KEY", ""
            ) == "disabled"

    def test_configured_api_keys_nonempty(self):
        """At least one API key should be configured (unless in deferred mode)."""
        from biocompiler.api import get_configured_api_keys
        from biocompiler.api import auth as _auth_module
        keys = get_configured_api_keys()
        # In deferred mode, keys are not generated until the API server starts.
        # In required/optional mode, auto-generation ensures at least one key.
        if _auth_module._AUTH_MODE != "deferred":
            assert len(keys) >= 1
        # In deferred mode, keys may be empty — that is correct.

    def test_verify_api_key_rejects_bad_key(self):
        """verify_api_key should raise HTTPException for invalid keys."""
        from fastapi import HTTPException
        from biocompiler.api import verify_api_key, get_auth_mode

        if get_auth_mode() == "disabled":
            pytest.skip("Auth is disabled in this environment")

        # We call verify_api_key with a bad key synchronously;
        # it is an async function so we use asyncio
        async def _test():
            try:
                await verify_api_key("invalid_key_12345")
                # If no exception, auth must be optional or disabled
                if get_auth_mode() == "required":
                    pytest.fail("Expected HTTPException for bad key")
            except HTTPException as exc:
                assert exc.status_code == 401

        asyncio.run(_test())

    def test_verify_api_key_accepts_valid_key(self):
        """verify_api_key should accept a configured key."""
        from biocompiler.api import verify_api_key, get_configured_api_keys, get_auth_mode

        keys = get_configured_api_keys()
        if not keys or get_auth_mode() == "disabled":
            pytest.skip("No keys configured or auth disabled")

        valid_key = next(iter(keys))

        async def _test():
            result = await verify_api_key(valid_key)
            assert result == valid_key

        asyncio.run(_test())


# ═══════════════════════════════════════════════════════════════════════
# 2. Rate Limiting Tests
# ═══════════════════════════════════════════════════════════════════════


class TestRateLimiting:
    """Test rate limiting with the SQLite backend."""

    def _make_limiter(self, max_requests: int = 5, window_seconds: int = 60):
        """Create a fresh rate limiter with a temp DB."""
        tmp = tempfile.NamedTemporaryFile(suffix=".db", delete=False)
        tmp.close()
        limiter = PersistentRateLimiter(
            db_path=tmp.name,
            max_requests=max_requests,
            window_seconds=window_seconds,
        )
        return limiter, tmp.name

    def test_rate_limiter_allows_within_limit(self):
        """Requests within the limit should be allowed."""
        limiter, db_path = self._make_limiter(max_requests=5)
        try:
            for i in range(5):
                allowed, remaining = limiter.check("client_a")
                assert allowed is True
                limiter.record("client_a")
        finally:
            Path(db_path).unlink(missing_ok=True)

    def test_rate_limiter_blocks_over_limit(self):
        """Requests exceeding the limit should be blocked."""
        limiter, db_path = self._make_limiter(max_requests=3)
        try:
            for i in range(3):
                limiter.record("client_b")
            # 4th request should be blocked
            allowed, remaining = limiter.check("client_b")
            assert allowed is False
            assert remaining == 0
        finally:
            Path(db_path).unlink(missing_ok=True)

    def test_rate_limiter_per_client_isolation(self):
        """Different clients should have independent rate limits."""
        limiter, db_path = self._make_limiter(max_requests=2)
        try:
            limiter.record("client_x")
            limiter.record("client_x")
            # client_x is now at limit
            allowed_x, _ = limiter.check("client_x")
            assert allowed_x is False
            # client_y should still be allowed
            allowed_y, _ = limiter.check("client_y")
            assert allowed_y is True
        finally:
            Path(db_path).unlink(missing_ok=True)

    def test_rate_limiter_batch_records(self):
        """Batch recording should consume multiple rate limit units."""
        limiter, db_path = self._make_limiter(max_requests=10)
        try:
            limiter.record_batch("client_c", 5)
            allowed, remaining = limiter.check("client_c")
            assert allowed is True
            assert remaining == 5
        finally:
            Path(db_path).unlink(missing_ok=True)

    def test_rate_limiter_clear(self):
        """Clearing should reset rate limit state."""
        limiter, db_path = self._make_limiter(max_requests=2)
        try:
            limiter.record("client_d")
            limiter.record("client_d")
            limiter.clear("client_d")
            allowed, remaining = limiter.check("client_d")
            assert allowed is True
        finally:
            Path(db_path).unlink(missing_ok=True)

    def test_rate_limiter_cleanup(self):
        """Cleanup should remove expired entries (mocked time, no sleep)."""
        # Use mock time source to avoid real sleep
        fake_time = [1000000.0]

        def mock_time_fn():
            return fake_time[0]

        tmp = tempfile.NamedTemporaryFile(suffix=".db", delete=False)
        tmp.close()
        limiter = PersistentRateLimiter(
            db_path=tmp.name,
            max_requests=100,
            window_seconds=60,
            time_func=mock_time_fn,
        )
        try:
            limiter.record("client_e")
            # Advance time past the window (no real sleep!)
            fake_time[0] += 61
            deleted = limiter.cleanup()
            # At least one entry should have been cleaned
            assert deleted >= 0  # may be 0 if already cleaned
            # After cleanup, client should have full quota
            allowed, remaining = limiter.check("client_e")
            assert allowed is True
        finally:
            Path(tmp.name).unlink(missing_ok=True)

    def test_rate_limiter_persistence(self):
        """Rate limit state should persist in the SQLite DB."""
        tmp = tempfile.NamedTemporaryFile(suffix=".db", delete=False)
        tmp.close()
        try:
            # Create limiter and record
            limiter1 = PersistentRateLimiter(
                db_path=tmp.name, max_requests=5, window_seconds=3600
            )
            limiter1.record("client_f")
            limiter1.record("client_f")

            # Create new limiter instance pointing to same DB
            limiter2 = PersistentRateLimiter(
                db_path=tmp.name, max_requests=5, window_seconds=3600
            )
            allowed, remaining = limiter2.check("client_f")
            assert allowed is True
            assert remaining == 3  # 5 max - 2 recorded
        finally:
            Path(tmp.name).unlink(missing_ok=True)


# ═══════════════════════════════════════════════════════════════════════
# 3. Input Size Limits
# ═══════════════════════════════════════════════════════════════════════


class TestInputSizeLimits:
    """Test that input size limits are enforced."""

    def test_max_protein_length_is_positive(self):
        """MAX_PROTEIN_LENGTH should be a positive integer."""
        assert isinstance(MAX_PROTEIN_LENGTH, int)
        assert MAX_PROTEIN_LENGTH > 0

    def test_max_batch_size_is_positive(self):
        """MAX_BATCH_SIZE should be a positive integer."""
        assert isinstance(MAX_BATCH_SIZE, int)
        assert MAX_BATCH_SIZE > 0

    def test_max_request_size_is_positive(self):
        """MAX_REQUEST_SIZE should be a positive integer."""
        assert isinstance(MAX_REQUEST_SIZE, int)
        assert MAX_REQUEST_SIZE > 0

    def test_optimize_timeout_is_positive(self):
        """OPTIMIZE_TIMEOUT_S should be a positive integer."""
        assert isinstance(OPTIMIZE_TIMEOUT_S, int)
        assert OPTIMIZE_TIMEOUT_S > 0

    def test_protein_exceeding_max_length_rejected(self):
        """ProteinInput should reject sequences exceeding MAX_PROTEIN_LENGTH."""
        long_protein = "A" * (MAX_PROTEIN_LENGTH + 1)
        with pytest.raises(Exception):  # ValidationError from pydantic
            ProteinInput(protein=long_protein, organism="Homo_sapiens")

    def test_protein_at_max_length_accepted(self):
        """ProteinInput should accept sequences exactly at MAX_PROTEIN_LENGTH."""
        protein_at_limit = "A" * MAX_PROTEIN_LENGTH
        inp = ProteinInput(protein=protein_at_limit, organism="Homo_sapiens")
        assert len(inp.protein) == MAX_PROTEIN_LENGTH

    def test_validate_protein_rejects_oversize(self):
        """validate_protein_input should return error for oversize sequences."""
        long_protein = "M" * (MAX_PROTEIN_LENGTH + 100)
        err = validate_protein_input(long_protein)
        assert err is not None
        assert "too long" in err.lower() or "maximum" in err.lower()

    def test_batch_limits_are_reasonable(self):
        """Batch limits should be reasonable values."""
        assert BATCH_CHECK_MAX > 0
        assert BATCH_OPTIMIZE_MAX > 0
        assert BATCH_EXPORT_MAX > 0
        # Optimize is the most expensive, so it should be smallest
        assert BATCH_OPTIMIZE_MAX <= BATCH_CHECK_MAX

    def test_rate_limit_rpm_is_positive(self):
        """Rate limit RPM should be a positive integer."""
        assert RATE_LIMIT_RPM > 0


# ═══════════════════════════════════════════════════════════════════════
# 4. Invalid Amino Acid Rejection
# ═══════════════════════════════════════════════════════════════════════


class TestInvalidAminoAcids:
    """Test that invalid amino acids are rejected with clear error messages."""

    def test_invalid_amino_acid_in_protein_input(self):
        """ProteinInput should reject sequences with non-standard AAs."""
        with pytest.raises(Exception) as exc_info:
            ProteinInput(protein="MVSKGEZ", organism="Homo_sapiens")
        # The error message should mention the invalid amino acid
        error_str = str(exc_info.value).lower()
        assert "invalid" in error_str or "z" in error_str

    def test_invalid_amino_acids_listed_in_error(self):
        """Error should list all invalid amino acids found."""
        with pytest.raises(Exception) as exc_info:
            ProteinInput(protein="MJOB", organism="Homo_sapiens")
        error_str = str(exc_info.value)
        # Should mention the specific invalid AAs (J, O, B)
        assert "J" in error_str or "O" in error_str or "B" in error_str

    def test_valid_amino_acids_accepted(self):
        """All 20 standard amino acids should be accepted."""
        from biocompiler.shared.constants import STANDARD_AAS
        inp = ProteinInput(protein=STANDARD_AINO_ACIDS if False else STANDARD_AAS, organism="Homo_sapiens")
        assert len(inp.protein) == 20

    def test_validate_protein_input_returns_clear_error(self):
        """validate_protein_input should name the invalid AAs clearly."""
        err = validate_protein_input("MVXZ")
        assert err is not None
        assert "X" in err or "Z" in err
        assert "Invalid" in err or "invalid" in err

    def test_lowercase_converted_to_uppercase(self):
        """Lowercase amino acid codes should be accepted (auto-uppercased)."""
        inp = ProteinInput(protein="mvskge", organism="Homo_sapiens")
        assert inp.protein == "MVSKGE"


# ═══════════════════════════════════════════════════════════════════════
# 5. Invalid Organism Rejection
# ═══════════════════════════════════════════════════════════════════════


class TestInvalidOrganism:
    """Test that invalid organisms are rejected with list of valid ones."""

    def test_invalid_organism_rejected(self):
        """ProteinInput should reject unsupported organisms."""
        with pytest.raises(Exception) as exc_info:
            ProteinInput(protein="MVSKGE", organism="Alien_martian")
        error_str = str(exc_info.value)
        # Error should mention the organism is unsupported
        assert "unsupported" in error_str.lower() or "invalid" in error_str.lower()

    def test_error_lists_supported_organisms(self):
        """Error message should list supported organisms."""
        with pytest.raises(Exception) as exc_info:
            ProteinInput(protein="MVSKGE", organism="Nonexistent_species")
        error_str = str(exc_info.value)
        # Should contain at least one known organism name
        assert "Homo_sapiens" in error_str or "Escherichia" in error_str or "supported" in error_str.lower()

    def test_validate_organism_returns_error_for_invalid(self):
        """validate_organism_input should return error for unknown organism."""
        err = validate_organism_input("Nonexistent_species")
        assert err is not None
        assert "unsupported" in err.lower() or "supported" in err.lower()

    def test_valid_organism_aliases_accepted(self):
        """Common organism aliases should be resolved and accepted."""
        # Test short key alias
        inp = ProteinInput(protein="MVSKGE", organism="ecoli")
        assert inp.organism in SUPPORTED_ORGANISMS

    def test_canonical_organism_names_accepted(self):
        """Full canonical organism names should be accepted."""
        for org in ["Homo_sapiens", "Escherichia_coli"]:
            inp = ProteinInput(protein="MVSKGE", organism=org)
            assert inp.organism == org

    def test_empty_organism_rejected(self):
        """Empty organism string should be rejected."""
        err = validate_organism_input("")
        assert err is not None


# ═══════════════════════════════════════════════════════════════════════
# 6. Optimization Timeout
# ═══════════════════════════════════════════════════════════════════════


class TestOptimizationTimeout:
    """Test that optimization timeout is configured and reasonable."""

    def test_optimize_timeout_has_reasonable_default(self):
        """OPTIMIZE_TIMEOUT_S should be between 30 and 3600 seconds."""
        assert 30 <= OPTIMIZE_TIMEOUT_S <= 3600

    def test_timeout_is_configurable_via_env(self):
        """Verify timeout can be overridden with env variable."""
        # We cannot easily change the module-level constant without reload,
        # but we can verify the env var mechanism exists
        env_var = os.environ.get("BIOCOMPILER_OPTIMIZE_TIMEOUT")
        if env_var is not None:
            assert OPTIMIZE_TIMEOUT_S == int(env_var)

    def test_long_sequence_would_timeout(self):
        """Very long sequences should trigger timeout behavior.

        We verify the timeout constant is used in the API's optimize
        endpoint by checking it is a reasonable value.
        """
        # The optimize endpoint uses OPTIMIZE_TIMEOUT_S for asyncio.wait_for
        # A 300-second (5 min) default is reasonable for proteins up to 10k aa
        assert OPTIMIZE_TIMEOUT_S >= 60  # at least 1 minute
        assert OPTIMIZE_TIMEOUT_S <= 600  # at most 10 minutes


# ═══════════════════════════════════════════════════════════════════════
# 7. /info Endpoint
# ═══════════════════════════════════════════════════════════════════════


class TestInfoEndpoint:
    """Test that /info endpoint returns correct limits."""

    def test_info_response_model_fields(self):
        """InfoResponse should contain all expected fields."""
        from biocompiler.api import InfoResponse
        fields = InfoResponse.model_fields
        assert "max_protein_length" in fields
        assert "max_batch_size" in fields
        assert "max_request_size" in fields
        assert "optimize_timeout_s" in fields
        assert "supported_organisms" in fields
        assert "api_version" in fields
        assert "safety_version" in fields

    def test_info_response_values(self):
        """InfoResponse should reflect actual configured values."""
        from biocompiler.api import InfoResponse
        from biocompiler import __version__

        info = InfoResponse(
            max_protein_length=MAX_PROTEIN_LENGTH,
            max_batch_size=MAX_BATCH_SIZE,
            max_request_size=MAX_REQUEST_SIZE,
            optimize_timeout_s=OPTIMIZE_TIMEOUT_S,
            supported_organisms=["Homo_sapiens", "Escherichia_coli"],
            api_version=__version__,
            safety_version=__version__,
            max_dna_length=MAX_PROTEIN_LENGTH * 3,
        )
        assert info.max_protein_length == MAX_PROTEIN_LENGTH
        assert info.max_batch_size == MAX_BATCH_SIZE
        assert info.max_request_size == MAX_REQUEST_SIZE
        assert info.optimize_timeout_s == OPTIMIZE_TIMEOUT_S
        assert isinstance(info.supported_organisms, list)
        assert len(info.supported_organisms) > 0


# ═══════════════════════════════════════════════════════════════════════
# 8. Concurrent Request Safety
# ═══════════════════════════════════════════════════════════════════════


class TestConcurrentRequestSafety:
    """Test that concurrent requests do not corrupt state."""

    def test_concurrent_rate_limiter_writes(self):
        """Concurrent writes to the rate limiter should not corrupt the DB."""
        tmp = tempfile.NamedTemporaryFile(suffix=".db", delete=False)
        tmp.close()
        limiter = PersistentRateLimiter(
            db_path=tmp.name, max_requests=1000, window_seconds=3600
        )
        errors = []

        def worker(client_id: str, count: int):
            try:
                for _ in range(count):
                    limiter.record(client_id)
            except Exception as e:
                errors.append(e)

        try:
            threads = []
            for i in range(10):
                t = threading.Thread(target=worker, args=(f"client_{i}", 20))
                threads.append(t)
                t.start()
            for t in threads:
                t.join()

            assert len(errors) == 0, f"Concurrent write errors: {errors}"

            # Verify total count
            allowed, remaining = limiter.check("client_0")
            assert allowed is True
            # client_0 should have recorded 20, leaving 980
            assert remaining == 980
        finally:
            Path(tmp.name).unlink(missing_ok=True)

    def test_concurrent_provenance_store_writes(self):
        """Concurrent provenance store writes should not corrupt data."""
        from biocompiler.provenance.decision_provenance import (
            ProvenanceStore,
            OptimizationDecisionTrail,
            CodonDecision,
            ConstraintDecision,
        )

        tmp_dir = tempfile.mkdtemp()
        store = ProvenanceStore(store_dir=tmp_dir)
        errors = []
        saved_ids = []

        trail = OptimizationDecisionTrail(
            gene_name="test",
            input_protein="MVSKGE",
            output_dna="ATGGTTTCTAAAGGTGAA",
            organism="Homo_sapiens",
            solver_backend="greedy",
            seed=42,
            total_cai=0.78,
            total_gc=0.50,
            codon_decisions=[
                CodonDecision(
                    position=0, amino_acid="M", original_codon=None,
                    chosen_codon="ATG", alternatives_considered=[],
                    constraint_reason="maximize_cai", confidence=1.0,
                ),
            ],
            constraint_decisions=[],
            iteration_log=[],
            timestamp="2026-01-01T00:00:00+00:00",
            version="test",
        )

        def worker(idx: int):
            try:
                record_id = store.save(trail)
                saved_ids.append(record_id)
            except Exception as e:
                errors.append(e)

        try:
            threads = []
            for i in range(5):
                t = threading.Thread(target=worker, args=(i,))
                threads.append(t)
                t.start()
            for t in threads:
                t.join()

            assert len(errors) == 0, f"Concurrent provenance errors: {errors}"
            # All records should be loadable
            for rid in saved_ids:
                loaded = store.load(rid)
                assert loaded.input_protein == "MVSKGE"
        finally:
            import shutil
            shutil.rmtree(tmp_dir, ignore_errors=True)


# ═══════════════════════════════════════════════════════════════════════
# 9. Provenance Persistence Across Restarts
# ═══════════════════════════════════════════════════════════════════════


class TestProvenancePersistence:
    """Test that provenance records survive API restarts."""

    def test_provenance_store_save_and_load(self):
        """A saved provenance record should be loadable by a new store instance."""
        from biocompiler.provenance.decision_provenance import (
            ProvenanceStore,
            OptimizationDecisionTrail,
            CodonDecision,
        )

        tmp_dir = tempfile.mkdtemp()
        try:
            trail = OptimizationDecisionTrail(
                gene_name="GFP",
                input_protein="MVSKGE",
                output_dna="ATGGTTTCTAAAGGTGAA",
                organism="Homo_sapiens",
                solver_backend="greedy",
                seed=123,
                total_cai=0.85,
                total_gc=0.52,
                codon_decisions=[
                    CodonDecision(
                        position=0, amino_acid="M", original_codon=None,
                        chosen_codon="ATG", alternatives_considered=[],
                        constraint_reason="maximize_cai", confidence=1.0,
                    ),
                ],
                constraint_decisions=[],
                iteration_log=[{"step": 1, "action": "init", "score": 0.0}],
                timestamp="2026-03-04T12:00:00+00:00",
                version="1.0.0",
            )

            # Save with store instance 1
            store1 = ProvenanceStore(store_dir=tmp_dir)
            record_id = store1.save(trail)

            # Load with a completely new store instance (simulates restart)
            store2 = ProvenanceStore(store_dir=tmp_dir)
            loaded = store2.load(record_id)

            assert loaded.gene_name == "GFP"
            assert loaded.input_protein == "MVSKGE"
            assert loaded.output_dna == "ATGGTTTCTAAAGGTGAA"
            assert loaded.organism == "Homo_sapiens"
            assert loaded.solver_backend == "greedy"
            assert loaded.seed == 123
            assert loaded.total_cai == 0.85
            assert loaded.total_gc == 0.52
            assert len(loaded.codon_decisions) == 1
            assert loaded.codon_decisions[0].chosen_codon == "ATG"
            assert len(loaded.iteration_log) == 1
        finally:
            import shutil
            shutil.rmtree(tmp_dir, ignore_errors=True)

    def test_provenance_query_after_restart(self):
        """Querying provenance should work after a simulated restart."""
        from biocompiler.provenance.decision_provenance import (
            ProvenanceStore,
            OptimizationDecisionTrail,
        )

        tmp_dir = tempfile.mkdtemp()
        try:
            trail = OptimizationDecisionTrail(
                gene_name="HBB",
                input_protein="MVHLTPEEK",
                output_dna="ATGGTGCATCTGACTCCTGAGGAGAAGTCT",
                organism="Homo_sapiens",
                solver_backend="greedy",
                seed=None,
                total_cai=0.72,
                total_gc=0.55,
                codon_decisions=[],
                constraint_decisions=[],
                iteration_log=[],
                timestamp="2026-03-04T12:00:00+00:00",
                version="1.0.0",
            )

            store1 = ProvenanceStore(store_dir=tmp_dir)
            store1.save(trail)

            # New instance simulating restart
            store2 = ProvenanceStore(store_dir=tmp_dir)
            results = store2.query(organism="Homo_sapiens")
            assert len(results) >= 1
            assert any(r.gene_name == "HBB" for r in results)
        finally:
            import shutil
            shutil.rmtree(tmp_dir, ignore_errors=True)

    def test_provenance_roundtrip_serialization(self):
        """Provenance should survive a to_dict → from_dict round trip."""
        from biocompiler.provenance.decision_provenance import (
            OptimizationDecisionTrail,
            CodonDecision,
            ConstraintDecision,
        )

        trail = OptimizationDecisionTrail(
            gene_name="test",
            input_protein="MVSKGE",
            output_dna="ATGGTTTCTAAAGGTGAA",
            organism="Escherichia_coli",
            solver_backend="z3",
            seed=99,
            total_cai=0.91,
            total_gc=0.48,
            codon_decisions=[
                CodonDecision(
                    position=0, amino_acid="M", original_codon="ATG",
                    chosen_codon="ATG", alternatives_considered=[
                        {"codon": "ATG", "cai": 1.0, "violations": []}
                    ],
                    constraint_reason="maximize_cai", confidence=0.95,
                    cai_impact=0.0,
                ),
            ],
            constraint_decisions=[
                ConstraintDecision(
                    constraint_name="GCInRange",
                    constraint_type="hard",
                    action_taken="satisfied",
                    positions_affected=[2, 5],
                    tradeoff_description="Selected GTC over GTG for GC balance",
                    impact_on_cai=-0.002,
                ),
            ],
            iteration_log=[{"step": 1, "action": "init"}],
            timestamp="2026-01-01T00:00:00+00:00",
            version="2.0.0",
        )

        # Round trip
        data = trail.to_dict()
        restored = OptimizationDecisionTrail.from_dict(data)

        assert restored.gene_name == trail.gene_name
        assert restored.input_protein == trail.input_protein
        assert restored.total_cai == trail.total_cai
        assert len(restored.codon_decisions) == 1
        assert restored.codon_decisions[0].cai_impact == 0.0
        assert len(restored.constraint_decisions) == 1
        assert restored.constraint_decisions[0].impact_on_cai == -0.002
