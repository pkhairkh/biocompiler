"""Integration tests for MHC adapter fallback chain (Task F4.4).

Tests that the mhcflurry_adapter.py works correctly in offline mode by
using the precomputed databases as fallback when MHCflurry models are
not installed.

Fallback chain:
  1. MHCflurry (neural network) — confidence 1.0
  2. NetMHCpan (web API) — confidence 1.0
  3. Precomputed database — confidence 0.7
  4. PSSM fallback — confidence 0.5

Every test is designed to run in a fully offline environment where
neither MHCflurry nor NetMHCpan are installed.
"""
from __future__ import annotations

import logging
import sys
import unittest
from dataclasses import replace
from unittest.mock import MagicMock, patch

# Ensure src is on the path for direct imports
sys.path.insert(0, "src")

from biocompiler.immunogenicity import (
    MHCBindingResult,
    classify_binding,
    score_peptide_pssm,
    binding_score_to_ic50,
)
from biocompiler.mhcflurry_adapter import (
    CONFIDENCE_MHCFLURRY,
    CONFIDENCE_NETMHCPAN,
    CONFIDENCE_PRECOMPUTED,
    CONFIDENCE_PSSM,
    MHCflurryClient,
    clear_cache,
    is_mhcflurry_available,
    is_netmhcpan_available,
    predict_binding,
    predict_batch,
)


class TestMHCBindingResultFields(unittest.TestCase):
    """Verify MHCBindingResult has the new rank and confidence fields."""

    def test_default_rank_is_none(self):
        result = MHCBindingResult(
            allele="HLA-A*02:01",
            peptide="LLFGYPVYV",
            start_position=0,
            end_position=8,
            binding_score=0.8,
            ic50_nm=50.0,
            binding_class="strong_binder",
            anchor_residues={},
            anchor_scores={},
        )
        self.assertIsNone(result.rank)

    def test_default_confidence_is_pssm(self):
        result = MHCBindingResult(
            allele="HLA-A*02:01",
            peptide="LLFGYPVYV",
            start_position=0,
            end_position=8,
            binding_score=0.5,
            ic50_nm=500.0,
            binding_class="moderate_binder",
            anchor_residues={},
            anchor_scores={},
        )
        self.assertEqual(result.confidence, 0.5)

    def test_explicit_rank_and_confidence(self):
        result = MHCBindingResult(
            allele="HLA-A*02:01",
            peptide="LLFGYPVYV",
            start_position=0,
            end_position=8,
            binding_score=0.9,
            ic50_nm=20.0,
            binding_class="strong_binder",
            anchor_residues={},
            anchor_scores={},
            method="mhcflurry",
            rank=0.5,
            confidence=1.0,
        )
        self.assertEqual(result.rank, 0.5)
        self.assertEqual(result.confidence, 1.0)

    def test_backward_compatible_without_new_fields(self):
        """Old code that doesn't pass rank/confidence should still work."""
        result = MHCBindingResult(
            allele="HLA-A*02:01",
            peptide="GILGFVFTL",
            start_position=0,
            end_position=8,
            binding_score=0.6,
            ic50_nm=300.0,
            binding_class="moderate_binder",
            anchor_residues={1: "I", 8: "L"},
            anchor_scores={1: 1.4, 8: 1.8},
            method="pssm",
        )
        self.assertIsNone(result.rank)
        self.assertEqual(result.confidence, 0.5)


class TestOfflinePredictBinding(unittest.TestCase):
    """Test predict_binding() works when MHCflurry is NOT installed."""

    def setUp(self):
        clear_cache()

    @patch("biocompiler.mhcflurry_adapter.is_mhcflurry_available", return_value=False)
    @patch("biocompiler.mhcflurry_adapter.is_netmhcpan_available", return_value=False)
    def test_predict_binding_returns_result_offline(self, mock_netmhcpan, mock_mhcflurry):
        """When neither MHCflurry nor NetMHCpan is available, we get a result."""
        result = predict_binding("HLA-A*02:01", "LLFGYPVYV")
        self.assertIsInstance(result, MHCBindingResult)
        self.assertEqual(result.allele, "HLA-A*02:01")
        self.assertEqual(result.peptide, "LLFGYPVYV")

    @patch("biocompiler.mhcflurry_adapter.is_mhcflurry_available", return_value=False)
    @patch("biocompiler.mhcflurry_adapter.is_netmhcpan_available", return_value=False)
    def test_predict_binding_uses_pssm_fallback(self, mock_netmhcpan, mock_mhcflurry):
        """In offline mode, PSSM fallback is used for known alleles."""
        result = predict_binding("HLA-A*02:01", "LLFGYPVYV")
        self.assertIn(result.method, ("pssm_fallback", "precomputed_lookup"))
        self.assertLessEqual(result.confidence, CONFIDENCE_PRECOMPUTED)

    @patch("biocompiler.mhcflurry_adapter.is_mhcflurry_available", return_value=False)
    @patch("biocompiler.mhcflurry_adapter.is_netmhcpan_available", return_value=False)
    def test_predict_binding_unknown_allele(self, mock_netmhcpan, mock_mhcflurry):
        """Unknown alleles should still return a result, not crash."""
        result = predict_binding("HLA-A*99:99", "LLFGYPVYV")
        self.assertIsInstance(result, MHCBindingResult)
        self.assertEqual(result.method, "pssm_fallback")
        self.assertEqual(result.confidence, CONFIDENCE_PSSM)

    @patch("biocompiler.mhcflurry_adapter.is_mhcflurry_available", return_value=False)
    @patch("biocompiler.mhcflurry_adapter.is_netmhcpan_available", return_value=False)
    def test_predict_binding_invalid_peptide_raises(self, mock_netmhcpan, mock_mhcflurry):
        """Invalid peptides should raise ValueError, not crash silently."""
        with self.assertRaises(ValueError):
            predict_binding("HLA-A*02:01", "123BADSEQ")


class TestPrecomputedDatabaseFallback(unittest.TestCase):
    """Test that precomputed database is used as fallback."""

    def setUp(self):
        clear_cache()

    @patch("biocompiler.mhcflurry_adapter.is_mhcflurry_available", return_value=False)
    @patch("biocompiler.mhcflurry_adapter.is_netmhcpan_available", return_value=False)
    def test_precomputed_lookup_for_known_epitope(self, mock_netmhcpan, mock_mhcflurry):
        """Known epitopes in the precomputed DB should be found via precomputed_lookup."""
        client = MHCflurryClient(allow_offline_fallback=True)

        # Mock the precomputed DB to return a result
        mock_entry = MagicMock()
        mock_entry.peptide = "LLFGYPVYV"
        mock_entry.binding_score = 1.0
        mock_entry.ic50_nm = 28.12
        mock_entry.binding_class = "strong_binder"
        mock_entry.anchor_residues = {0: "L", 1: "L", 8: "V"}
        mock_entry.anchor_scores = {0: 1.2, 1: 2.0, 8: 2.0}
        mock_entry.source = "known_epitope"

        with patch.object(client, "_lookup_precomputed", return_value=mock_entry):
            result = client.predict_binding("LLFGYPVYV", "HLA-A*02:01")
            self.assertEqual(result.method, "precomputed_lookup")
            self.assertEqual(result.confidence, CONFIDENCE_PRECOMPUTED)
            self.assertAlmostEqual(result.ic50_nm, 28.12, places=1)

    @patch("biocompiler.mhcflurry_adapter.is_mhcflurry_available", return_value=False)
    @patch("biocompiler.mhcflurry_adapter.is_netmhcpan_available", return_value=False)
    def test_precomputed_miss_falls_to_pssm(self, mock_netmhcpan, mock_mhcflurry):
        """When precomputed DB has no entry, PSSM fallback is used."""
        client = MHCflurryClient(allow_offline_fallback=True)

        with patch.object(client, "_lookup_precomputed", return_value=None):
            result = client.predict_binding("LLFGYPVYV", "HLA-A*02:01")
            self.assertEqual(result.method, "pssm_fallback")
            self.assertEqual(result.confidence, CONFIDENCE_PSSM)


class TestPSSMFallback(unittest.TestCase):
    """Test that PSSM fallback works when precomputed data is not available."""

    def setUp(self):
        clear_cache()

    @patch("biocompiler.mhcflurry_adapter.is_mhcflurry_available", return_value=False)
    @patch("biocompiler.mhcflurry_adapter.is_netmhcpan_available", return_value=False)
    def test_pssm_fallback_for_known_allele(self, mock_netmhcpan, mock_mhcflurry):
        """PSSM fallback should produce reasonable results for alleles with PSSMs."""
        client = MHCflurryClient(allow_offline_fallback=True)

        with patch.object(client, "_lookup_precomputed", return_value=None):
            result = client.predict_binding("LLFGYPVYV", "HLA-A*02:01")
            self.assertEqual(result.method, "pssm_fallback")
            self.assertEqual(result.confidence, CONFIDENCE_PSSM)
            # LLFGYPVYV has L at pos 1 (score 2.0) and V at pos 8 (score 2.0)
            # which are strong anchors — should bind
            self.assertIsNotNone(result.ic50_nm)

    @patch("biocompiler.mhcflurry_adapter.is_mhcflurry_available", return_value=False)
    @patch("biocompiler.mhcflurry_adapter.is_netmhcpan_available", return_value=False)
    def test_pssm_fallback_for_unknown_allele(self, mock_netmhcpan, mock_mhcflurry):
        """Unknown alleles should return non-binder from PSSM fallback."""
        client = MHCflurryClient(allow_offline_fallback=True)

        with patch.object(client, "_lookup_precomputed", return_value=None):
            result = client.predict_binding("LLFGYPVYV", "HLA-XYZ*99:99")
            self.assertEqual(result.method, "pssm_fallback")
            self.assertEqual(result.confidence, CONFIDENCE_PSSM)
            # No PSSM for unknown allele — ic50 should be very high (non-binder)
            self.assertGreater(result.ic50_nm, 5000.0)

    @patch("biocompiler.mhcflurry_adapter.is_mhcflurry_available", return_value=False)
    @patch("biocompiler.mhcflurry_adapter.is_netmhcpan_available", return_value=False)
    def test_pssm_fallback_never_crashes(self, mock_netmhcpan, mock_mhcflurry):
        """PSSM fallback should never crash, even with edge-case inputs."""
        client = MHCflurryClient(allow_offline_fallback=True)

        with patch.object(client, "_lookup_precomputed", return_value=None):
            # Short peptide (8-mer, different from 9-mer PSSM)
            result = client.predict_binding("SIINFEKL", "HLA-A*02:01")
            self.assertIsInstance(result, MHCBindingResult)
            # 8-mer doesn't match 9-mer PSSM → should return non-binder
            self.assertIn(result.method, ("pssm_fallback",))
            self.assertEqual(result.confidence, CONFIDENCE_PSSM)


class TestBatchPrediction(unittest.TestCase):
    """Test batch prediction via predict_batch()."""

    def setUp(self):
        clear_cache()

    @patch("biocompiler.mhcflurry_adapter.is_mhcflurry_available", return_value=False)
    @patch("biocompiler.mhcflurry_adapter.is_netmhcpan_available", return_value=False)
    def test_predict_batch_returns_one_per_peptide(self, mock_netmhcpan, mock_mhcflurry):
        """predict_batch should return one result per input peptide."""
        peptides = ["LLFGYPVYV", "SIINFEKL", "GILGFVFTL"]
        results = predict_batch("HLA-A*02:01", peptides)
        self.assertEqual(len(results), len(peptides))
        for r in results:
            self.assertIsInstance(r, MHCBindingResult)

    @patch("biocompiler.mhcflurry_adapter.is_mhcflurry_available", return_value=False)
    @patch("biocompiler.mhcflurry_adapter.is_netmhcpan_available", return_value=False)
    def test_predict_batch_preserves_peptide_order(self, mock_netmhcpan, mock_mhcflurry):
        """Results should be in the same order as input peptides."""
        peptides = ["GILGFVFTL", "LLFGYPVYV"]
        results = predict_batch("HLA-A*02:01", peptides)
        self.assertEqual(results[0].peptide, "GILGFVFTL")
        self.assertEqual(results[1].peptide, "LLFGYPVYV")

    @patch("biocompiler.mhcflurry_adapter.is_mhcflurry_available", return_value=False)
    @patch("biocompiler.mhcflurry_adapter.is_netmhcpan_available", return_value=False)
    def test_predict_batch_empty_list(self, mock_netmhcpan, mock_mhcflurry):
        """Empty peptide list should return empty results."""
        results = predict_batch("HLA-A*02:01", [])
        self.assertEqual(len(results), 0)


class TestConfidenceLevels(unittest.TestCase):
    """Test that confidence levels are set correctly for each fallback tier."""

    def setUp(self):
        clear_cache()

    def test_confidence_constants(self):
        """Verify the expected confidence constant values."""
        self.assertEqual(CONFIDENCE_MHCFLURRY, 1.0)
        self.assertEqual(CONFIDENCE_NETMHCPAN, 1.0)
        self.assertEqual(CONFIDENCE_PRECOMPUTED, 0.7)
        self.assertEqual(CONFIDENCE_PSSM, 0.5)

    @patch("biocompiler.mhcflurry_adapter.is_mhcflurry_available", return_value=False)
    @patch("biocompiler.mhcflurry_adapter.is_netmhcpan_available", return_value=False)
    def test_precomputed_confidence(self, mock_netmhcpan, mock_mhcflurry):
        """Precomputed lookup should set confidence to 0.7."""
        client = MHCflurryClient(allow_offline_fallback=True)

        mock_entry = MagicMock()
        mock_entry.peptide = "TESTPEPTI"
        mock_entry.binding_score = 0.8
        mock_entry.ic50_nm = 100.0
        mock_entry.binding_class = "moderate_binder"
        mock_entry.anchor_residues = {}
        mock_entry.anchor_scores = {}
        mock_entry.source = "pssm_predicted"

        with patch.object(client, "_lookup_precomputed", return_value=mock_entry):
            result = client.predict_binding("TESTPEPTI", "HLA-A*02:01")
            self.assertEqual(result.method, "precomputed_lookup")
            self.assertEqual(result.confidence, CONFIDENCE_PRECOMPUTED)

    @patch("biocompiler.mhcflurry_adapter.is_mhcflurry_available", return_value=False)
    @patch("biocompiler.mhcflurry_adapter.is_netmhcpan_available", return_value=False)
    def test_pssm_confidence(self, mock_netmhcpan, mock_mhcflurry):
        """PSSM fallback should set confidence to 0.5."""
        client = MHCflurryClient(allow_offline_fallback=True)

        with patch.object(client, "_lookup_precomputed", return_value=None):
            result = client.predict_binding("LLFGYPVYV", "HLA-A*02:01")
            self.assertEqual(result.method, "pssm_fallback")
            self.assertEqual(result.confidence, CONFIDENCE_PSSM)

    def test_mhcflurry_result_has_high_confidence(self):
        """MHCflurry results created via _mhcflurry_result_to_binding_result
        should have confidence 1.0."""
        from biocompiler.mhcflurry_adapter import _mhcflurry_result_to_binding_result

        result = _mhcflurry_result_to_binding_result(
            peptide="LLFGYPVYV",
            allele="HLA-A*02:01",
            start_position=0,
            end_position=8,
            ic50_nm=50.0,
            method="mhcflurry",
        )
        self.assertEqual(result.confidence, CONFIDENCE_MHCFLURRY)

    def test_netmhcpan_result_confidence(self):
        """NetMHCpan results should have confidence 1.0."""
        result = MHCBindingResult(
            allele="HLA-A*02:01",
            peptide="LLFGYPVYV",
            start_position=0,
            end_position=8,
            binding_score=0.9,
            ic50_nm=30.0,
            binding_class="strong_binder",
            anchor_residues={},
            anchor_scores={},
            method="netmhcpan",
            confidence=CONFIDENCE_NETMHCPAN,
        )
        self.assertEqual(result.confidence, 1.0)


class TestNeverCrashes(unittest.TestCase):
    """Verify the adapter NEVER crashes due to missing models."""

    def setUp(self):
        clear_cache()

    @patch("biocompiler.mhcflurry_adapter.is_mhcflurry_available", return_value=False)
    @patch("biocompiler.mhcflurry_adapter.is_netmhcpan_available", return_value=False)
    def test_predict_binding_no_mhcflurry(self, mock_netmhcpan, mock_mhcflurry):
        """predict_binding should not crash when MHCflurry is missing."""
        result = predict_binding("HLA-A*02:01", "LLFGYPVYV")
        self.assertIsInstance(result, MHCBindingResult)

    def test_client_load_models_does_not_raise(self):
        """_load_models should return False instead of raising on ImportError."""
        client = MHCflurryClient(allow_offline_fallback=True)
        # _load_models should never raise — it returns bool
        result = client._load_models()
        self.assertIsInstance(result, bool)

    @patch("biocompiler.mhcflurry_adapter.is_mhcflurry_available", return_value=False)
    @patch("biocompiler.mhcflurry_adapter.is_netmhcpan_available", return_value=False)
    def test_client_predict_binding_no_crash(self, mock_netmhcpan, mock_mhcflurry):
        """MHCflurryClient.predict_binding should never crash in offline mode."""
        client = MHCflurryClient(allow_offline_fallback=True)
        result = client.predict_binding("LLFGYPVYV", "HLA-A*02:01")
        self.assertIsInstance(result, MHCBindingResult)

    @patch("biocompiler.mhcflurry_adapter.is_mhcflurry_available", return_value=False)
    @patch("biocompiler.mhcflurry_adapter.is_netmhcpan_available", return_value=False)
    def test_client_batch_predict_no_crash(self, mock_netmhcpan, mock_mhcflurry):
        """MHCflurryClient.batch_predict should never crash in offline mode."""
        client = MHCflurryClient(allow_offline_fallback=True)
        results = client.batch_predict(
            "MAGRSGDLDAIIRYVKQLR",
            alleles=["HLA-A*02:01"],
            epitope_lengths=[9],
        )
        self.assertIsInstance(results, list)
        # Should have results for each 9-mer in the protein
        self.assertGreater(len(results), 0)
        for r in results:
            self.assertIsInstance(r, MHCBindingResult)


class TestFallbackLogging(unittest.TestCase):
    """Verify that fallback transitions are properly logged."""

    def setUp(self):
        clear_cache()

    @patch("biocompiler.mhcflurry_adapter.is_mhcflurry_available", return_value=False)
    @patch("biocompiler.mhcflurry_adapter.is_netmhcpan_available", return_value=False)
    def test_mhcflurry_unavailable_logs_info(self, mock_netmhcpan, mock_mhcflurry):
        """When MHCflurry is unavailable, an INFO log should be emitted."""
        client = MHCflurryClient(allow_offline_fallback=True)

        with self.assertLogs("biocompiler.mhcflurry_adapter", level="INFO") as cm:
            with patch.object(client, "_lookup_precomputed", return_value=None):
                client.predict_binding("LLFGYPVYV", "HLA-A*02:01")

        # Should log about fallback transitions
        log_messages = "\n".join(cm.output)
        # Should mention either MHCflurry unavailable or precomputed/PSSM fallback
        self.assertTrue(
            any("MHCflurry" in msg for msg in cm.output),
            f"Expected MHCflurry fallback log, got: {cm.output}",
        )

    @patch("biocompiler.mhcflurry_adapter.is_mhcflurry_available", return_value=False)
    @patch("biocompiler.mhcflurry_adapter.is_netmhcpan_available", return_value=False)
    def test_precomputed_miss_logs_pssm_fallback(self, mock_netmhcpan, mock_mhcflurry):
        """When precomputed DB has no match, PSSM fallback should be logged."""
        client = MHCflurryClient(allow_offline_fallback=True)

        with self.assertLogs("biocompiler.mhcflurry_adapter", level="INFO") as cm:
            with patch.object(client, "_lookup_precomputed", return_value=None):
                client.predict_binding("LLFGYPVYV", "HLA-A*02:01")

        log_messages = "\n".join(cm.output)
        self.assertTrue(
            any("PSSM" in msg or "precomputed" in msg for msg in cm.output),
            f"Expected PSSM fallback log, got: {cm.output}",
        )


class TestModuleLevelFunctions(unittest.TestCase):
    """Test module-level predict_binding and predict_batch convenience functions."""

    def setUp(self):
        clear_cache()

    def test_predict_binding_is_callable(self):
        """Module-level predict_binding should be a callable function."""
        self.assertTrue(callable(predict_binding))

    def test_predict_batch_is_callable(self):
        """Module-level predict_batch should be a callable function."""
        self.assertTrue(callable(predict_batch))

    @patch("biocompiler.mhcflurry_adapter.is_mhcflurry_available", return_value=False)
    @patch("biocompiler.mhcflurry_adapter.is_netmhcpan_available", return_value=False)
    def test_predict_binding_returns_mhc_result(self, mock_netmhcpan, mock_mhcflurry):
        """predict_binding should return MHCBindingResult."""
        result = predict_binding("HLA-A*02:01", "LLFGYPVYV")
        self.assertIsInstance(result, MHCBindingResult)
        self.assertEqual(result.allele, "HLA-A*02:01")

    @patch("biocompiler.mhcflurry_adapter.is_mhcflurry_available", return_value=False)
    @patch("biocompiler.mhcflurry_adapter.is_netmhcpan_available", return_value=False)
    def test_predict_batch_returns_list(self, mock_netmhcpan, mock_mhcflurry):
        """predict_batch should return a list of MHCBindingResult."""
        results = predict_batch("HLA-A*02:01", ["LLFGYPVYV", "GILGFVFTL"])
        self.assertIsInstance(results, list)
        self.assertEqual(len(results), 2)
        for r in results:
            self.assertIsInstance(r, MHCBindingResult)


class TestNetMHCpanFallback(unittest.TestCase):
    """Test that NetMHCpan is tried as second fallback."""

    def setUp(self):
        clear_cache()

    @patch("biocompiler.mhcflurry_adapter.is_mhcflurry_available", return_value=False)
    @patch("biocompiler.mhcflurry_adapter.is_netmhcpan_available", return_value=True)
    def test_netmhcpan_tried_when_mhcflurry_unavailable(self, mock_netmhcpan_avail, mock_mhcflurry_avail):
        """When MHCflurry is unavailable but NetMHCpan is, try NetMHCpan."""
        client = MHCflurryClient(allow_offline_fallback=True)

        # Mock the NetMHCpan prediction to return None (not available after all)
        with patch.object(client, "_try_netmhcpan_prediction", return_value=None) as mock_netmhcpan_pred:
            with patch.object(client, "_lookup_precomputed", return_value=None):
                result = client.predict_binding("LLFGYPVYV", "HLA-A*02:01")
                # Should have tried NetMHCpan
                mock_netmhcpan_pred.assert_called_once_with("LLFGYPVYV", "HLA-A*02:01")
                # Falls through to PSSM
                self.assertEqual(result.method, "pssm_fallback")

    @patch("biocompiler.mhcflurry_adapter.is_mhcflurry_available", return_value=False)
    @patch("biocompiler.mhcflurry_adapter.is_netmhcpan_available", return_value=True)
    def test_netmhcpan_success_uses_netmhcpan_method(self, mock_netmhcpan_avail, mock_mhcflurry_avail):
        """When NetMHCpan succeeds, the result should have method='netmhcpan'."""
        client = MHCflurryClient(allow_offline_fallback=True)

        netmhcpan_result = MHCBindingResult(
            allele="HLA-A*02:01",
            peptide="LLFGYPVYV",
            start_position=0,
            end_position=8,
            binding_score=0.85,
            ic50_nm=35.0,
            binding_class="strong_binder",
            anchor_residues={},
            anchor_scores={},
            method="netmhcpan",
            confidence=CONFIDENCE_NETMHCPAN,
        )
        with patch.object(client, "_try_netmhcpan_prediction", return_value=netmhcpan_result):
            result = client.predict_binding("LLFGYPVYV", "HLA-A*02:01")
            self.assertEqual(result.method, "netmhcpan")
            self.assertEqual(result.confidence, CONFIDENCE_NETMHCPAN)


if __name__ == "__main__":
    unittest.main()
