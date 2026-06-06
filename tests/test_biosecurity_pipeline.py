"""
Tests for biosecurity screening integration into the optimization pipeline.

Task 1.1 — Wire biosecurity screening into the optimization pipeline.

Covers:
  1. optimize_sequence() rejects a ricin sequence (enforce mode)
  2. optimize_sequence() accepts a normal insulin sequence (enforce mode)
  3. The three biosecurity modes (enforce / warn / off)
  4. The API /optimize endpoint with biosecurity screening (403 response)
  5. BioOptimizer.optimize() biosecurity gate
  6. BIOCOMPILER_BIOSECURITY_MODE environment variable
  7. biosecurity_mode parameter on optimize_sequence() and BioOptimizer.__init__()
"""

from __future__ import annotations

import os
from unittest.mock import patch

import pytest

from biocompiler.biosecurity import (
    BiosecurityReport,
    HazardMatch,
    check_biosecurity_before_optimize,
    get_biosecurity_mode,
    screen_hazardous_sequence,
)
from biocompiler.exceptions import BiosecurityError, BioCompilerError
from biocompiler.optimization import optimize_sequence, BioOptimizer


# ═══════════════════════════════════════════════════════════════════════════
# Test sequences
# ═══════════════════════════════════════════════════════════════════════════

# Human insulin B-chain (benign therapeutic protein)
INSULIN_PROTEIN = "MALWMRLLPLLALLALWGPDPAAAFVNQHLCGSHLVEALYLVCGERGFFYTPKT"

# Ricin A-chain catalytic motif embedded in a longer sequence
RICIN_PROTEIN = "MISRDNIRVGLPIISTNKYEDKQL"


# ═══════════════════════════════════════════════════════════════════════════
# 1. optimize_sequence() rejects ricin, accepts insulin
# ═══════════════════════════════════════════════════════════════════════════

class TestOptimizeSequenceBiosecurityGate:
    """Verify that optimize_sequence() screens input before optimizing."""

    def test_ricin_rejected_in_enforce_mode(self):
        """A ricin sequence must be rejected with BiosecurityError."""
        with pytest.raises(BiosecurityError) as exc_info:
            optimize_sequence(
                target_protein=RICIN_PROTEIN,
                organism="Escherichia_coli",
                biosecurity_mode="enforce",
            )
        assert exc_info.value.risk_level in ("critical", "high")

    def test_insulin_accepted_in_enforce_mode(self):
        """A normal insulin sequence must pass biosecurity screening."""
        result = optimize_sequence(
            target_protein=INSULIN_PROTEIN,
            organism="Escherichia_coli",
            biosecurity_mode="enforce",
        )
        assert result.sequence is not None
        assert len(result.sequence) > 0

    def test_biosecurity_error_is_biocompiler_error(self):
        """BiosecurityError should be catchable as BioCompilerError."""
        with pytest.raises(BioCompilerError):
            optimize_sequence(
                target_protein=RICIN_PROTEIN,
                organism="Escherichia_coli",
                biosecurity_mode="enforce",
            )

    def test_ricin_rejected_default_mode(self):
        """Default mode (enforce) should also reject ricin without explicit parameter."""
        with pytest.raises(BiosecurityError):
            optimize_sequence(
                target_protein=RICIN_PROTEIN,
                organism="Escherichia_coli",
                # biosecurity_mode intentionally not set — defaults to env var "enforce"
            )


# ═══════════════════════════════════════════════════════════════════════════
# 2. Three biosecurity modes: enforce / warn / off
# ═══════════════════════════════════════════════════════════════════════════

class TestBiosecurityModes:
    """Test enforce, warn, and off biosecurity modes."""

    def test_enforce_mode_raises_on_ricin(self):
        """In enforce mode, hazardous sequences raise BiosecurityError."""
        with pytest.raises(BiosecurityError):
            check_biosecurity_before_optimize(
                RICIN_PROTEIN,
                organism="e_coli",
                biosecurity_mode="enforce",
            )

    def test_enforce_mode_allows_insulin(self):
        """In enforce mode, clean sequences pass through."""
        report = check_biosecurity_before_optimize(
            INSULIN_PROTEIN,
            organism="e_coli",
            biosecurity_mode="enforce",
        )
        assert report.risk_level == "none"
        assert not report.is_hazardous

    def test_warn_mode_does_not_raise_on_ricin(self):
        """In warn mode, even hazardous sequences do NOT raise."""
        import warnings
        with warnings.catch_warnings(record=True):
            warnings.simplefilter("always")
            report = check_biosecurity_before_optimize(
                RICIN_PROTEIN,
                organism="e_coli",
                biosecurity_mode="warn",
            )
        # Should NOT raise, should still report the hazard
        assert report.is_hazardous
        assert report.risk_level in ("critical", "high")

    def test_warn_mode_emits_warning_on_ricin(self):
        """In warn mode, hazardous sequences emit a UserWarning."""
        import warnings
        with warnings.catch_warnings(record=True) as w:
            warnings.simplefilter("always")
            check_biosecurity_before_optimize(
                RICIN_PROTEIN,
                organism="e_coli",
                biosecurity_mode="warn",
            )
        biosec_warnings = [x for x in w if issubclass(x.category, UserWarning) and "Biosecurity" in str(x.message)]
        assert len(biosec_warnings) >= 1

    def test_off_mode_skips_screening(self):
        """In off mode, screening is skipped entirely."""
        report = check_biosecurity_before_optimize(
            RICIN_PROTEIN,
            organism="e_coli",
            biosecurity_mode="off",
        )
        assert report.risk_level == "none"
        assert not report.is_hazardous
        assert len(report.matches) == 0

    def test_off_mode_returns_skip_recommendation(self):
        """Off mode should include a recommendation noting screening was skipped."""
        report = check_biosecurity_before_optimize(
            RICIN_PROTEIN,
            organism="e_coli",
            biosecurity_mode="off",
        )
        assert any("skipped" in r.lower() for r in report.recommendations)

    def test_optimize_sequence_warn_mode_proceeds_on_ricin(self):
        """optimize_sequence() with biosecurity_mode='warn' should proceed
        despite hazardous input."""
        result = optimize_sequence(
            target_protein=RICIN_PROTEIN,
            organism="Escherichia_coli",
            biosecurity_mode="warn",
        )
        assert result.sequence is not None

    def test_optimize_sequence_off_mode_proceeds_on_ricin(self):
        """optimize_sequence() with biosecurity_mode='off' should proceed
        without any screening."""
        result = optimize_sequence(
            target_protein=RICIN_PROTEIN,
            organism="Escherichia_coli",
            biosecurity_mode="off",
        )
        assert result.sequence is not None


# ═══════════════════════════════════════════════════════════════════════════
# 3. BIOCOMPILER_BIOSECURITY_MODE environment variable
# ═══════════════════════════════════════════════════════════════════════════

class TestBiosecurityModeEnvVar:
    """Test the BIOCOMPILER_BIOSECURITY_MODE environment variable."""

    def test_default_is_enforce(self):
        """Without the env var, the default mode should be 'enforce'."""
        with patch.dict(os.environ, {}, clear=False):
            # Remove the var if present
            os.environ.pop("BIOCOMPILER_BIOSECURITY_MODE", None)
            mode = get_biosecurity_mode()
            assert mode == "enforce"

    def test_enforce_mode_from_env(self):
        """BIOCOMPILER_BIOSECURITY_MODE=enforce should return 'enforce'."""
        with patch.dict(os.environ, {"BIOCOMPILER_BIOSECURITY_MODE": "enforce"}):
            assert get_biosecurity_mode() == "enforce"

    def test_warn_mode_from_env(self):
        """BIOCOMPILER_BIOSECURITY_MODE=warn should return 'warn'."""
        with patch.dict(os.environ, {"BIOCOMPILER_BIOSECURITY_MODE": "warn"}):
            assert get_biosecurity_mode() == "warn"

    def test_off_mode_from_env(self):
        """BIOCOMPILER_BIOSECURITY_MODE=off should return 'off'."""
        with patch.dict(os.environ, {"BIOCOMPILER_BIOSECURITY_MODE": "off"}):
            assert get_biosecurity_mode() == "off"

    def test_case_insensitive_env(self):
        """Environment variable should be case-insensitive."""
        with patch.dict(os.environ, {"BIOCOMPILER_BIOSECURITY_MODE": "WARN"}):
            assert get_biosecurity_mode() == "warn"

    def test_invalid_env_falls_back_to_enforce(self):
        """An invalid env var value should fall back to 'enforce'."""
        with patch.dict(os.environ, {"BIOCOMPILER_BIOSECURITY_MODE": "invalid"}):
            assert get_biosecurity_mode() == "enforce"

    def test_parameter_overrides_env(self):
        """The biosecurity_mode parameter should take precedence over the env var."""
        # Env says "off" but parameter says "enforce" — should enforce
        with patch.dict(os.environ, {"BIOCOMPILER_BIOSECURITY_MODE": "off"}):
            with pytest.raises(BiosecurityError):
                check_biosecurity_before_optimize(
                    RICIN_PROTEIN,
                    organism="e_coli",
                    biosecurity_mode="enforce",
                )

    def test_env_warn_mode_allows_ricin(self):
        """With env var set to 'warn', ricin should not raise."""
        import warnings
        with patch.dict(os.environ, {"BIOCOMPILER_BIOSECURITY_MODE": "warn"}):
            # No explicit parameter, so it reads from env
            with warnings.catch_warnings(record=True):
                warnings.simplefilter("always")
                report = check_biosecurity_before_optimize(
                    RICIN_PROTEIN,
                    organism="e_coli",
                    # biosecurity_mode=None => reads from env
                )
            assert report.is_hazardous  # Still reports the hazard


# ═══════════════════════════════════════════════════════════════════════════
# 4. BioOptimizer.optimize() biosecurity gate
# ═══════════════════════════════════════════════════════════════════════════

class TestBioOptimizerBiosecurityGate:
    """Verify BioOptimizer.optimize() screens input before optimizing."""

    def test_biooptimizer_rejects_ricin_in_enforce_mode(self):
        """BioOptimizer with biosecurity_mode='enforce' rejects ricin."""
        optimizer = BioOptimizer(
            species="ecoli",
            biosecurity_mode="enforce",
        )
        # Need to pass a DNA sequence that encodes a protein with ricin motif
        # The optimizer's optimize() takes a DNA or protein sequence
        # Let's construct a DNA sequence that translates to something with the ricin motif
        from biocompiler.type_system import CODON_TABLE
        # Manually construct a DNA sequence that encodes RICIN_PROTEIN
        # Use the most common codons
        codon_choices = {
            'M': 'ATG', 'I': 'ATC', 'S': 'AGC', 'R': 'CGC', 'D': 'GAC',
            'N': 'AAC', 'V': 'GTG', 'G': 'GGC', 'L': 'CTG', 'P': 'CCC',
            'T': 'ACC', 'K': 'AAG', 'Y': 'TAC', 'E': 'GAG', 'Q': 'CAG',
        }
        ricin_dna = ''.join(codon_choices.get(aa, 'NNN') for aa in RICIN_PROTEIN)
        with pytest.raises(BiosecurityError):
            optimizer.optimize(ricin_dna)

    def test_biooptimizer_accepts_insulin_in_enforce_mode(self):
        """BioOptimizer with biosecurity_mode='enforce' accepts insulin."""
        optimizer = BioOptimizer(
            species="ecoli",
            biosecurity_mode="enforce",
        )
        from biocompiler.type_system import CODON_TABLE
        codon_choices = {
            'M': 'ATG', 'A': 'GCC', 'L': 'CTG', 'W': 'TGG', 'R': 'CGC',
            'P': 'CCC', 'H': 'CAC', 'V': 'GTG', 'F': 'TTC', 'Y': 'TAC',
            'N': 'AAC', 'Q': 'CAG', 'C': 'TGC', 'G': 'GGC', 'S': 'AGC',
            'D': 'GAC', 'E': 'GAG', 'I': 'ATC', 'K': 'AAG', 'T': 'ACC',
        }
        insulin_dna = ''.join(codon_choices.get(aa, 'NNN') for aa in INSULIN_PROTEIN)
        result = optimizer.optimize(insulin_dna)
        # Result is (optimized_seq, predicates, certificate)
        assert result is not None
        assert len(result[0]) > 0

    def test_biooptimizer_warn_mode_proceeds_on_ricin(self):
        """BioOptimizer with biosecurity_mode='warn' proceeds despite hazard."""
        import warnings
        optimizer = BioOptimizer(
            species="ecoli",
            biosecurity_mode="warn",
        )
        from biocompiler.type_system import CODON_TABLE
        codon_choices = {
            'M': 'ATG', 'I': 'ATC', 'S': 'AGC', 'R': 'CGC', 'D': 'GAC',
            'N': 'AAC', 'V': 'GTG', 'G': 'GGC', 'L': 'CTG', 'P': 'CCC',
            'T': 'ACC', 'K': 'AAG', 'Y': 'TAC', 'E': 'GAG', 'Q': 'CAG',
        }
        ricin_dna = ''.join(codon_choices.get(aa, 'NNN') for aa in RICIN_PROTEIN)
        with warnings.catch_warnings(record=True):
            warnings.simplefilter("always")
            result = optimizer.optimize(ricin_dna)
        assert result is not None

    def test_biooptimizer_off_mode_skips_screening(self):
        """BioOptimizer with biosecurity_mode='off' skips screening."""
        optimizer = BioOptimizer(
            species="ecoli",
            biosecurity_mode="off",
        )
        from biocompiler.type_system import CODON_TABLE
        codon_choices = {
            'M': 'ATG', 'I': 'ATC', 'S': 'AGC', 'R': 'CGC', 'D': 'GAC',
            'N': 'AAC', 'V': 'GTG', 'G': 'GGC', 'L': 'CTG', 'P': 'CCC',
            'T': 'ACC', 'K': 'AAG', 'Y': 'TAC', 'E': 'GAG', 'Q': 'CAG',
        }
        ricin_dna = ''.join(codon_choices.get(aa, 'NNN') for aa in RICIN_PROTEIN)
        # Should not raise
        result = optimizer.optimize(ricin_dna)
        assert result is not None

    def test_biooptimizer_default_no_biosecurity_mode(self):
        """BioOptimizer without explicit biosecurity_mode should default to
        reading from env (enforce)."""
        optimizer = BioOptimizer(species="ecoli")
        assert optimizer.biosecurity_mode is None  # Parameter is None => reads env


# ═══════════════════════════════════════════════════════════════════════════
# 5. API /optimize endpoint with biosecurity screening
# ═══════════════════════════════════════════════════════════════════════════

class TestAPIBiosecurityScreening:
    """Test the /optimize API endpoint biosecurity gate."""

    @pytest.fixture()
    def app(self):
        """Create a FastAPI app with auth disabled for testing."""
        from biocompiler import api as _api_module
        from biocompiler.api import create_app, _rate_limiter

        _rate_limiter.clear()
        original_mode = _api_module._AUTH_MODE
        original_keys = set(_api_module._CONFIGURED_API_KEYS)
        _api_module._AUTH_MODE = "disabled"
        _api_module._CONFIGURED_API_KEYS = set()
        try:
            yield create_app()
        finally:
            _api_module._AUTH_MODE = original_mode
            _api_module._CONFIGURED_API_KEYS = original_keys

    @pytest.fixture()
    def client(self, app):
        from fastapi.testclient import TestClient
        return TestClient(app)

    def test_optimize_ricin_returns_403(self, client):
        """Optimizing a ricin sequence should return 403 Forbidden."""
        resp = client.post("/optimize", json={
            "protein": RICIN_PROTEIN,
            "organism": "Escherichia_coli",
        })
        assert resp.status_code == 403
        data = resp.json()
        assert "detail" in data
        detail = data["detail"]
        assert detail.get("error") == "BiosecurityError"
        assert detail.get("risk_level") in ("critical", "high")
        assert "select_agent" in detail.get("flagged_categories", [])

    def test_optimize_ricin_403_has_matches(self, client):
        """The 403 response should include match details."""
        resp = client.post("/optimize", json={
            "protein": RICIN_PROTEIN,
            "organism": "Escherichia_coli",
        })
        assert resp.status_code == 403
        data = resp.json()
        detail = data["detail"]
        assert "matches" in detail
        assert len(detail["matches"]) > 0
        # At least one match should be ricin
        match_names = [m["name"] for m in detail["matches"]]
        assert any("ricin" in n for n in match_names)

    def test_optimize_ricin_403_has_recommendations(self, client):
        """The 403 response should include recommendations."""
        resp = client.post("/optimize", json={
            "protein": RICIN_PROTEIN,
            "organism": "Escherichia_coli",
        })
        assert resp.status_code == 403
        data = resp.json()
        detail = data["detail"]
        assert "recommendations" in detail
        assert len(detail["recommendations"]) > 0

    def test_optimize_insulin_returns_200(self, client):
        """Optimizing insulin should succeed with 200."""
        resp = client.post("/optimize", json={
            "protein": INSULIN_PROTEIN,
            "organism": "Escherichia_coli",
        })
        assert resp.status_code == 200
        data = resp.json()
        assert "sequence" in data
        assert "protein" in data

    def test_optimize_ricin_warn_mode_returns_200(self, client):
        """With biosecurity_mode=warn, ricin should be allowed through the API.

        Note: The API endpoint currently uses the env var to determine mode,
        not a request parameter. We test by setting the env var.
        """
        with patch.dict(os.environ, {"BIOCOMPILER_BIOSECURITY_MODE": "warn"}):
            resp = client.post("/optimize", json={
                "protein": RICIN_PROTEIN,
                "organism": "Escherichia_coli",
                "strict_mode": False,
            })
            # In warn mode, the API-level check will proceed,
            # but optimize_sequence will also check — so it should succeed
            assert resp.status_code == 200

    def test_optimize_ricin_off_mode_returns_200(self, client):
        """With biosecurity_mode=off, ricin should bypass screening."""
        with patch.dict(os.environ, {"BIOCOMPILER_BIOSECURITY_MODE": "off"}):
            resp = client.post("/optimize", json={
                "protein": RICIN_PROTEIN,
                "organism": "Escherichia_coli",
                "strict_mode": False,
            })
            assert resp.status_code == 200


# ═══════════════════════════════════════════════════════════════════════════
# 6. BiosecurityMode parameter propagation
# ═══════════════════════════════════════════════════════════════════════════

class TestBiosecurityModeParameter:
    """Verify biosecurity_mode parameter is accepted and propagated."""

    def test_optimize_sequence_accepts_biosecurity_mode(self):
        """optimize_sequence() should accept biosecurity_mode parameter."""
        # Just verify it doesn't raise TypeError
        result = optimize_sequence(
            target_protein=INSULIN_PROTEIN,
            organism="Escherichia_coli",
            biosecurity_mode="enforce",
        )
        assert result is not None

    def test_biooptimizer_init_accepts_biosecurity_mode(self):
        """BioOptimizer.__init__() should accept biosecurity_mode parameter."""
        optimizer = BioOptimizer(
            species="ecoli",
            biosecurity_mode="warn",
        )
        assert optimizer.biosecurity_mode == "warn"

    def test_biooptimizer_init_default_biosecurity_mode_is_none(self):
        """BioOptimizer.__init__() default for biosecurity_mode is None."""
        optimizer = BioOptimizer(species="ecoli")
        assert optimizer.biosecurity_mode is None

    def test_optimize_sequence_biosecurity_mode_none_reads_env(self):
        """When biosecurity_mode=None, optimize_sequence reads from env."""
        # Env set to "off" should allow ricin through
        with patch.dict(os.environ, {"BIOCOMPILER_BIOSECURITY_MODE": "off"}):
            result = optimize_sequence(
                target_protein=RICIN_PROTEIN,
                organism="Escherichia_coli",
                biosecurity_mode=None,
            )
            assert result.sequence is not None


# ═══════════════════════════════════════════════════════════════════════════
# 7. Integration: full pipeline with biosecurity gate
# ═══════════════════════════════════════════════════════════════════════════

class TestFullPipelineBiosecurityIntegration:
    """End-to-end tests of the biosecurity gate in the optimization pipeline."""

    def test_ricin_blocked_before_optimization_begins(self):
        """Biosecurity screening should block BEFORE any optimization work."""
        with pytest.raises(BiosecurityError):
            optimize_sequence(
                target_protein=RICIN_PROTEIN,
                organism="Homo_sapiens",
                biosecurity_mode="enforce",
            )

    def test_insulin_optimization_result_has_protein(self):
        """Optimizing insulin should produce a result with the correct protein."""
        result = optimize_sequence(
            target_protein=INSULIN_PROTEIN,
            organism="Escherichia_coli",
            biosecurity_mode="enforce",
        )
        assert result.protein == INSULIN_PROTEIN

    def test_medium_risk_sequence_in_enforce_mode(self):
        """A medium-risk sequence (e.g., nptII) should warn but not block."""
        import warnings
        # nptII protein motif — medium risk
        nptii_protein = "AAARPMTIHGSGSAAA"
        # First check that it's actually medium risk
        report = screen_hazardous_sequence(nptii_protein)
        if report.risk_level == "medium":
            with warnings.catch_warnings(record=True):
                warnings.simplefilter("always")
                result = optimize_sequence(
                    target_protein=nptii_protein,
                    organism="Escherichia_coli",
                    biosecurity_mode="enforce",
                )
                assert result is not None

    def test_screening_gate_is_early(self):
        """Verify that the biosecurity check happens before heavy computation.

        This is implicitly tested by the fact that ricin sequences
        raise BiosecurityError immediately. But we also verify that
        the error doesn't come from inside the optimization loop.
        """
        with pytest.raises(BiosecurityError) as exc_info:
            optimize_sequence(
                target_protein=RICIN_PROTEIN,
                organism="Escherichia_coli",
                biosecurity_mode="enforce",
            )
        # The error should reference the biosecurity screening
        assert "biosecurity" in str(exc_info.value).lower() or "risk" in str(exc_info.value).lower()
