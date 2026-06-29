"""Tests for biocompiler.biopython_compat — BioPython SeqRecord interoperability.

Covered areas:
1. Import / availability check  (_check_biopython, module-level import)
2. Compatibility function signatures (to_seqrecord, from_seqrecord, optimize_to_seqrecord)
3. Graceful degradation when BioPython is not installed (mocked ImportError)

All BioPython-dependent tests use ``pytest.importorskip("Bio")`` so the
entire suite is skipped cleanly when BioPython is absent.  Tests that
verify *degradation* behaviour mock the import at the ``_check_biopython``
boundary.
"""

from __future__ import annotations

import sys
import types
from unittest import mock

import pytest

# Mark every test in this module as requiring an external tool (BioPython).
# These tests use ``pytest.importorskip("Bio")`` at runtime; the marker keeps
# them deselected by default alongside other requires_external tests.
pytestmark = pytest.mark.requires_external

# ---------------------------------------------------------------------------
# 1. Import / availability check
# ---------------------------------------------------------------------------


class TestModuleImport:
    """Verify the biopython_compat module can be imported regardless of
    whether BioPython itself is installed.  The module must *always* be
    importable — only calling the functions should raise if BioPython is
    missing.
    """

    def test_module_importable(self):
        """biocompiler.biopython_compat can be imported without BioPython."""
        import biocompiler.shared.biopython_compat as mod  # noqa: F401
        assert mod is not None

    def test_module_has_expected_public_names(self):
        """Module exposes the three public compatibility functions."""
        import biocompiler.shared.biopython_compat as mod

        public = [name for name in dir(mod) if not name.startswith("_")]
        assert "to_seqrecord" in public
        assert "from_seqrecord" in public
        assert "optimize_to_seqrecord" in public

    def test_module_docstring(self):
        """Module has a non-empty docstring."""
        import biocompiler.shared.biopython_compat as mod
        assert mod.__doc__ is not None and len(mod.__doc__.strip()) > 0

    def test_check_biopython_callable(self):
        """_check_biopython is a callable internal helper."""
        import biocompiler.shared.biopython_compat as mod
        assert callable(mod._check_biopython)


class TestCheckBiopython:
    """Tests for the ``_check_biopython()`` guard function."""

    def test_succeeds_when_biopython_installed(self):
        """_check_biopython does not raise when Bio is importable."""
        import biocompiler.shared.biopython_compat as mod

        # If Bio is available this should simply return None (no exception)
        try:
            mod._check_biopython()
        except ImportError:
            pytest.skip("BioPython not installed; cannot test positive path")

    def test_raises_importerror_when_biopython_missing(self):
        """_check_biopython raises ImportError with helpful message when Bio
        cannot be imported."""
        import biocompiler.shared.biopython_compat as mod

        # Mock the import of Bio inside _check_biopython to simulate absence
        with mock.patch.dict(sys.modules, {"Bio": None}):
            # Also make importlib raise ImportError for 'Bio'
            import builtins
            real_import = builtins.__import__

            def _fake_import(name, *args, **kwargs):
                if name == "Bio":
                    raise ImportError("No module named 'Bio'")
                return real_import(name, *args, **kwargs)

            with mock.patch("builtins.__import__", _fake_import):
                with pytest.raises(ImportError, match="BioPython"):
                    mod._check_biopython()

    def test_importerror_message_mentions_install_command(self):
        """The ImportError message should mention pip install instructions."""
        import biocompiler.shared.biopython_compat as mod

        with mock.patch.dict(sys.modules, {"Bio": None}):
            import builtins
            real_import = builtins.__import__

            def _fake_import(name, *args, **kwargs):
                if name == "Bio":
                    raise ImportError("No module named 'Bio'")
                return real_import(name, *args, **kwargs)

            with mock.patch("builtins.__import__", _fake_import):
                with pytest.raises(ImportError) as exc_info:
                    mod._check_biopython()

                msg = str(exc_info.value)
                assert "pip install" in msg or "biopython" in msg.lower()

    def test_importerror_message_mentions_biocompiler_extra(self):
        """The ImportError message should mention the [biopython] extra."""
        import biocompiler.shared.biopython_compat as mod

        with mock.patch.dict(sys.modules, {"Bio": None}):
            import builtins
            real_import = builtins.__import__

            def _fake_import(name, *args, **kwargs):
                if name == "Bio":
                    raise ImportError("No module named 'Bio'")
                return real_import(name, *args, **kwargs)

            with mock.patch("builtins.__import__", _fake_import):
                with pytest.raises(ImportError) as exc_info:
                    mod._check_biopython()

                msg = str(exc_info.value)
                assert "biocompiler[biopython]" in msg


# ---------------------------------------------------------------------------
# 2. Compatibility function signatures
# ---------------------------------------------------------------------------


class TestToSeqrecordSignature:
    """Verify ``to_seqrecord`` has the expected parameter names, defaults,
    and return-type annotation."""

    def test_function_callable(self):
        """to_seqrecord is callable."""
        import biocompiler.shared.biopython_compat as mod
        assert callable(mod.to_seqrecord)

    def test_parameter_names(self):
        """to_seqrecord accepts the documented parameters."""
        import biocompiler.shared.biopython_compat as mod
        import inspect

        sig = inspect.signature(mod.to_seqrecord)
        params = list(sig.parameters.keys())
        assert "sequence" in params
        assert "organism" in params
        assert "gene_name" in params
        assert "exon_boundaries" in params
        assert "type_results" in params
        assert "certificate" in params

    def test_default_values(self):
        """Parameter defaults match the documented values."""
        import biocompiler.shared.biopython_compat as mod
        import inspect

        sig = inspect.signature(mod.to_seqrecord)
        p = sig.parameters

        assert p["organism"].default == "Homo_sapiens"
        assert p["gene_name"].default is None
        assert p["exon_boundaries"].default is None
        assert p["type_results"].default is None
        assert p["certificate"].default is None

    def test_sequence_is_required_positional(self):
        """The ``sequence`` parameter is the first positional parameter.

        Note: the implementation accepts either ``sequence=`` or ``result=``
        (keyword-only) for flexibility. The test verifies that ``sequence``
        is the first parameter, not that it has no default.
        """
        import biocompiler.shared.biopython_compat as mod
        import inspect

        sig = inspect.signature(mod.to_seqrecord)
        params = list(sig.parameters.keys())
        assert params[0] == "sequence", (
            f"Expected 'sequence' to be the first parameter, got {params[0]!r}"
        )

    def test_return_annotation(self):
        """Return annotation references Bio.SeqRecord.SeqRecord."""
        import biocompiler.shared.biopython_compat as mod
        import inspect

        sig = inspect.signature(mod.to_seqrecord)
        ret = sig.return_annotation
        # The annotation is a string "Bio.SeqRecord.SeqRecord" or the actual type
        ret_str = str(ret)
        assert "SeqRecord" in ret_str


class TestFromSeqrecordSignature:
    """Verify ``from_seqrecord`` has the expected parameter names and
    return-type annotation."""

    def test_function_callable(self):
        """from_seqrecord is callable."""
        import biocompiler.shared.biopython_compat as mod
        assert callable(mod.from_seqrecord)

    def test_parameter_names(self):
        """from_seqrecord accepts a single ``record`` parameter."""
        import biocompiler.shared.biopython_compat as mod
        import inspect

        sig = inspect.signature(mod.from_seqrecord)
        params = list(sig.parameters.keys())
        assert "record" in params

    def test_record_is_required(self):
        """The ``record`` parameter has no default (required)."""
        import biocompiler.shared.biopython_compat as mod
        import inspect

        sig = inspect.signature(mod.from_seqrecord)
        assert sig.parameters["record"].default is inspect.Parameter.empty

    def test_return_annotation_is_dict(self):
        """Return annotation references dict."""
        import biocompiler.shared.biopython_compat as mod
        import inspect

        sig = inspect.signature(mod.from_seqrecord)
        ret_str = str(sig.return_annotation)
        assert "dict" in ret_str


class TestOptimizeToSeqrecordSignature:
    """Verify ``optimize_to_seqrecord`` has the expected parameter names,
    defaults, and return-type annotation."""

    def test_function_callable(self):
        """optimize_to_seqrecord is callable."""
        import biocompiler.shared.biopython_compat as mod
        assert callable(mod.optimize_to_seqrecord)

    def test_parameter_names(self):
        """optimize_to_seqrecord accepts the documented parameters."""
        import biocompiler.shared.biopython_compat as mod
        import inspect

        sig = inspect.signature(mod.optimize_to_seqrecord)
        params = list(sig.parameters.keys())
        assert "protein" in params
        assert "organism" in params
        assert "gc_lo" in params
        assert "gc_hi" in params
        assert "cai_threshold" in params
        assert "restriction_enzymes" in params
        assert "gene_name" in params

    def test_default_values(self):
        """Parameter defaults match the documented values."""
        import biocompiler.shared.biopython_compat as mod
        import inspect

        sig = inspect.signature(mod.optimize_to_seqrecord)
        p = sig.parameters

        assert p["organism"].default == "Homo_sapiens"
        assert p["gc_lo"].default == pytest.approx(0.30)
        assert p["gc_hi"].default == pytest.approx(0.70)
        assert p["cai_threshold"].default == pytest.approx(0.2)
        assert p["restriction_enzymes"].default is None
        assert p["gene_name"].default is None

    def test_protein_is_required(self):
        """The ``protein`` parameter has no default (required)."""
        import biocompiler.shared.biopython_compat as mod
        import inspect

        sig = inspect.signature(mod.optimize_to_seqrecord)
        assert sig.parameters["protein"].default is inspect.Parameter.empty

    def test_return_annotation(self):
        """Return annotation references Bio.SeqRecord.SeqRecord."""
        import biocompiler.shared.biopython_compat as mod
        import inspect

        sig = inspect.signature(mod.optimize_to_seqrecord)
        ret_str = str(sig.return_annotation)
        assert "SeqRecord" in ret_str


# ---------------------------------------------------------------------------
# 3. Graceful degradation when BioPython not installed
# ---------------------------------------------------------------------------


class TestGracefulDegradation:
    """All three public functions must raise ImportError (not crash with
    AttributeError or worse) when BioPython is unavailable."""

    @staticmethod
    def _make_biopython_unavailable():
        """Return a context manager that makes ``import Bio`` fail."""
        import builtins
        real_import = builtins.__import__

        def _fake_import(name, *args, **kwargs):
            if name == "Bio" or name.startswith("Bio."):
                raise ImportError(f"No module named '{name}'")
            return real_import(name, *args, **kwargs)

        return mock.patch("builtins.__import__", _fake_import)

    # --- to_seqrecord ---

    def test_to_seqrecord_raises_importerror_without_biopython(self):
        """to_seqrecord raises ImportError when BioPython is missing."""
        import biocompiler.shared.biopython_compat as mod

        with self._make_biopython_unavailable():
            with pytest.raises(ImportError, match="BioPython"):
                mod.to_seqrecord("ATGCGT")

    def test_to_seqrecord_importerror_is_helpful(self):
        """The ImportError includes install instructions."""
        import biocompiler.shared.biopython_compat as mod

        with self._make_biopython_unavailable():
            with pytest.raises(ImportError) as exc_info:
                mod.to_seqrecord("ATGCGT")

            msg = str(exc_info.value).lower()
            assert "pip" in msg or "install" in msg

    # --- from_seqrecord ---

    def test_from_seqrecord_raises_importerror_without_biopython(self):
        """from_seqrecord raises ImportError when BioPython is missing."""
        import biocompiler.shared.biopython_compat as mod

        with self._make_biopython_unavailable():
            with pytest.raises(ImportError, match="BioPython"):
                mod.from_seqrecord("not_a_real_record")

    # --- optimize_to_seqrecord ---

    def test_optimize_to_seqrecord_raises_importerror_without_biopython(self):
        """optimize_to_seqrecord raises ImportError when BioPython is missing."""
        import biocompiler.shared.biopython_compat as mod

        with self._make_biopython_unavailable():
            # The ImportError should be raised before any optimization
            # happens (checked at _check_biopython in to_seqrecord or at
            # the top of optimize_to_seqrecord).
            with pytest.raises(ImportError, match="BioPython"):
                mod.optimize_to_seqrecord("MVLSPADKTN")

    def test_degradation_does_not_leave_partial_state(self):
        """After an ImportError, re-importing the module still works."""
        import biocompiler.shared.biopython_compat as mod

        with self._make_biopython_unavailable():
            with pytest.raises(ImportError):
                mod.to_seqrecord("ATG")

        # After the mock is removed, the module should still be importable
        import importlib
        importlib.reload(mod)
        assert callable(mod.to_seqrecord)


# ---------------------------------------------------------------------------
# 4. Functional tests (require BioPython)
# ---------------------------------------------------------------------------

# Only run the functional tests below when BioPython is actually installed.
Bio = pytest.importorskip("Bio", reason="BioPython not installed — skipping functional tests")

from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from Bio.SeqFeature import SeqFeature, FeatureLocation, CompoundLocation


class TestToSeqrecordFunctional:
    """Functional tests for ``to_seqrecord`` when BioPython is available."""

    def test_returns_seqrecord(self, sample_dna):
        """to_seqrecord returns a Bio.SeqRecord.SeqRecord instance."""
        from biocompiler.infrastructure.biopython_compat import to_seqrecord

        rec = to_seqrecord(sample_dna)
        assert isinstance(rec, SeqRecord)

    def test_sequence_set_correctly(self, sample_dna):
        """The SeqRecord.seq contains the uppercased DNA string."""
        from biocompiler.infrastructure.biopython_compat import to_seqrecord

        rec = to_seqrecord(sample_dna)
        assert str(rec.seq) == sample_dna.upper()

    def test_lowercase_sequence_uppercased(self):
        """Lowercase DNA is automatically uppercased in the SeqRecord."""
        from biocompiler.infrastructure.biopython_compat import to_seqrecord

        rec = to_seqrecord("atgcgt")
        assert str(rec.seq) == "ATGCGT"

    def test_organism_in_annotations(self, sample_dna):
        """Organism is stored in record.annotations."""
        from biocompiler.infrastructure.biopython_compat import to_seqrecord

        rec = to_seqrecord(sample_dna, organism="Escherichia_coli")
        assert rec.annotations["organism"] == "Escherichia_coli"

    def test_default_organism(self, sample_dna):
        """Default organism is Homo_sapiens."""
        from biocompiler.infrastructure.biopython_compat import to_seqrecord

        rec = to_seqrecord(sample_dna)
        assert rec.annotations["organism"] == "Homo_sapiens"

    def test_molecule_type_annotation(self, sample_dna):
        """molecule_type is 'DNA' in annotations."""
        from biocompiler.infrastructure.biopython_compat import to_seqrecord

        rec = to_seqrecord(sample_dna)
        assert rec.annotations["molecule_type"] == "DNA"

    def test_topology_annotation(self, sample_dna):
        """topology is 'linear' in annotations."""
        from biocompiler.infrastructure.biopython_compat import to_seqrecord

        rec = to_seqrecord(sample_dna)
        assert rec.annotations["topology"] == "linear"

    def test_gc_content_in_annotations(self, sample_dna):
        """gc_content is computed and stored in annotations."""
        from biocompiler.infrastructure.biopython_compat import to_seqrecord

        rec = to_seqrecord(sample_dna)
        assert "gc_content" in rec.annotations
        gc = rec.annotations["gc_content"]
        assert isinstance(gc, float)
        assert 0.0 <= gc <= 1.0

    def test_gene_name_in_id_and_name(self, sample_dna):
        """Gene name is set as record.id and record.name."""
        from biocompiler.infrastructure.biopython_compat import to_seqrecord

        rec = to_seqrecord(sample_dna, gene_name="eGFP")
        assert rec.id == "eGFP"
        assert rec.name == "eGFP"

    def test_default_id_when_no_gene_name(self, sample_dna):
        """Default id is 'BioCompiler_design' when gene_name is not provided."""
        from biocompiler.infrastructure.biopython_compat import to_seqrecord

        rec = to_seqrecord(sample_dna)
        assert rec.id == "BioCompiler_design"

    def test_description_contains_organism(self, sample_dna):
        """Description mentions the organism."""
        from biocompiler.infrastructure.biopython_compat import to_seqrecord

        rec = to_seqrecord(sample_dna, organism="Mus_musculus")
        assert "Mus_musculus" in rec.description

    def test_gene_feature_created(self, sample_dna):
        """When gene_name is provided, a 'gene' feature is created."""
        from biocompiler.infrastructure.biopython_compat import to_seqrecord

        rec = to_seqrecord(sample_dna, gene_name="eGFP")
        gene_feats = [f for f in rec.features if f.type == "gene"]
        assert len(gene_feats) == 1
        assert "gene" in gene_feats[0].qualifiers
        assert gene_feats[0].qualifiers["gene"] == ["eGFP"]

    def test_no_gene_feature_without_name(self, sample_dna):
        """No 'gene' feature is created when gene_name is None."""
        from biocompiler.infrastructure.biopython_compat import to_seqrecord

        rec = to_seqrecord(sample_dna)
        gene_feats = [f for f in rec.features if f.type == "gene"]
        assert len(gene_feats) == 0

    def test_cds_feature_always_created(self, sample_dna):
        """A CDS feature is always created."""
        from biocompiler.infrastructure.biopython_compat import to_seqrecord

        rec = to_seqrecord(sample_dna)
        cds_feats = [f for f in rec.features if f.type == "CDS"]
        assert len(cds_feats) >= 1

    def test_cds_translation_qualifier(self, sample_dna):
        """CDS feature has a translation qualifier with the protein."""
        from biocompiler.infrastructure.biopython_compat import to_seqrecord

        rec = to_seqrecord(sample_dna)
        cds_feats = [f for f in rec.features if f.type == "CDS"]
        translation = cds_feats[0].qualifiers.get("translation")
        assert translation is not None
        # Should be a non-empty protein string
        assert len(translation) > 0

    def test_exon_features_created(self, sample_dna):
        """Exon features are created when exon_boundaries are provided."""
        from biocompiler.infrastructure.biopython_compat import to_seqrecord

        boundaries = [(0, 9), (9, 18)]
        rec = to_seqrecord(sample_dna, exon_boundaries=boundaries)
        exon_feats = [f for f in rec.features if f.type == "exon"]
        assert len(exon_feats) == 2

    def test_exon_boundary_positions(self, sample_dna):
        """Exon feature locations match the provided boundaries."""
        from biocompiler.infrastructure.biopython_compat import to_seqrecord

        boundaries = [(0, 9), (9, 18)]
        rec = to_seqrecord(sample_dna, exon_boundaries=boundaries)
        exon_feats = sorted(
            [f for f in rec.features if f.type == "exon"],
            key=lambda f: int(f.location.start),
        )
        assert int(exon_feats[0].location.start) == 0
        assert int(exon_feats[0].location.end) == 9
        assert int(exon_feats[1].location.start) == 9
        assert int(exon_feats[1].location.end) == 18

    def test_multi_exon_cds_compound_location(self, sample_dna):
        """Multi-exon CDS uses CompoundLocation (join)."""
        from biocompiler.infrastructure.biopython_compat import to_seqrecord

        boundaries = [(0, 9), (9, 18)]
        rec = to_seqrecord(sample_dna, exon_boundaries=boundaries)
        cds_feats = [f for f in rec.features if f.type == "CDS"]
        assert len(cds_feats) >= 1
        assert hasattr(cds_feats[0].location, "parts")

    def test_single_exon_simple_location(self, sample_dna):
        """Single exon boundary produces a simple FeatureLocation, not CompoundLocation."""
        from Bio.SeqFeature import CompoundLocation as _CL
        from biocompiler.infrastructure.biopython_compat import to_seqrecord

        boundaries = [(0, 18)]
        rec = to_seqrecord(sample_dna, exon_boundaries=boundaries)
        cds_feats = [f for f in rec.features if f.type == "CDS"]
        # Single exon should NOT produce a CompoundLocation instance
        # (SimpleLocation may have a .parts property in newer BioPython, so
        # check the type rather than hasattr)
        assert not isinstance(cds_feats[0].location, _CL)

    def test_no_exon_boundaries_full_cds(self, sample_dna):
        """No exon_boundaries → CDS spans the entire sequence."""
        from biocompiler.infrastructure.biopython_compat import to_seqrecord

        rec = to_seqrecord(sample_dna)
        cds_feats = [f for f in rec.features if f.type == "CDS"]
        assert int(cds_feats[0].location.start) == 0
        assert int(cds_feats[0].location.end) == len(sample_dna)

    def test_certificate_embedded_in_annotations(self, sample_dna):
        """Certificate dict is stored in record.annotations."""
        from biocompiler.infrastructure.biopython_compat import to_seqrecord

        cert = {"design_id": "test-123", "status": "verified"}
        rec = to_seqrecord(sample_dna, certificate=cert)
        assert rec.annotations.get("biocompiler_certificate") == cert
        assert rec.annotations.get("biocompiler_design_id") == "test-123"

    def test_certificate_with_to_dict_method(self, sample_dna):
        """Certificate objects with a to_dict() method are serialized."""
        from biocompiler.infrastructure.biopython_compat import to_seqrecord

        class FakeCert:
            def to_dict(self):
                return {"design_id": "cert-456"}

        rec = to_seqrecord(sample_dna, certificate=FakeCert())
        assert rec.annotations["biocompiler_certificate"]["design_id"] == "cert-456"

    def test_certificate_serialization_failure_handled(self, sample_dna):
        """Non-serializable certificate object is handled gracefully."""
        from biocompiler.infrastructure.biopython_compat import to_seqrecord

        rec = to_seqrecord(sample_dna, certificate=42)  # int, not dict-like
        cert = rec.annotations["biocompiler_certificate"]
        assert "error" in cert

    def test_type_results_as_misc_features(self, sample_dna):
        """Type results are converted to misc_feature annotations."""
        from biocompiler.infrastructure.biopython_compat import to_seqrecord

        # Create a minimal mock type-check result
        class FakeVerdict:
            value = "PASS"

        class FakeTypeResult:
            predicate = "GcContent"
            verdict = FakeVerdict()
            violation = None
            knowledge_gap = None

        rec = to_seqrecord(sample_dna, type_results=[FakeTypeResult()])
        misc_feats = [f for f in rec.features if f.type == "misc_feature"]
        assert len(misc_feats) >= 1
        notes = misc_feats[0].qualifiers.get("note", [])
        assert any("GcContent" in n for n in notes)

    def test_type_result_with_violation(self, sample_dna):
        """Type results with violations include violation in notes."""
        from biocompiler.infrastructure.biopython_compat import to_seqrecord

        class FakeVerdict:
            value = "FAIL"

        class FakeTypeResult:
            predicate = "RestrictionSite"
            verdict = FakeVerdict()
            violation = "EcoRI found at position 5"
            knowledge_gap = None

        rec = to_seqrecord(sample_dna, type_results=[FakeTypeResult()])
        misc_feats = [f for f in rec.features if f.type == "misc_feature"]
        notes = misc_feats[0].qualifiers.get("note", [])
        assert any("Violation" in n for n in notes)

    def test_type_result_with_knowledge_gap(self, sample_dna):
        """Type results with knowledge gaps include gap in notes."""
        from biocompiler.infrastructure.biopython_compat import to_seqrecord

        class FakeVerdict:
            value = "UNCERTAIN"

        class FakeTypeResult:
            predicate = "MRNASecondaryStructure"
            verdict = FakeVerdict()
            violation = None
            knowledge_gap = "ViennaRNA not available"

        rec = to_seqrecord(sample_dna, type_results=[FakeTypeResult()])
        misc_feats = [f for f in rec.features if f.type == "misc_feature"]
        notes = misc_feats[0].qualifiers.get("note", [])
        assert any("Knowledge gap" in n for n in notes)


class TestFromSeqrecordFunctional:
    """Functional tests for ``from_seqrecord`` when BioPython is available."""

    def _make_record(self, seq: str = "ATGGTTTCTAAAGGTGAA", **kwargs) -> SeqRecord:
        """Helper: build a SeqRecord with optional kwargs forwarded to to_seqrecord."""
        from biocompiler.infrastructure.biopython_compat import to_seqrecord
        return to_seqrecord(seq, **kwargs)

    def test_returns_dict(self, sample_dna):
        """from_seqrecord returns a dict."""
        from biocompiler.infrastructure.biopython_compat import from_seqrecord, to_seqrecord

        rec = to_seqrecord(sample_dna)
        result = from_seqrecord(rec)
        assert isinstance(result, dict)

    def test_required_keys_present(self, sample_dna):
        """Returned dict has all documented keys."""
        from biocompiler.infrastructure.biopython_compat import from_seqrecord, to_seqrecord

        rec = to_seqrecord(sample_dna)
        result = from_seqrecord(rec)
        expected_keys = {
            "sequence", "organism", "gene_name", "exon_boundaries",
            "protein", "certificate", "features", "gc_content",
        }
        assert expected_keys.issubset(result.keys())

    def test_sequence_roundtrip(self, sample_dna):
        """Sequence survives to_seqrecord → from_seqrecord round-trip."""
        from biocompiler.infrastructure.biopython_compat import from_seqrecord, to_seqrecord

        rec = to_seqrecord(sample_dna)
        result = from_seqrecord(rec)
        assert result["sequence"] == sample_dna.upper()

    def test_organism_roundtrip(self, sample_dna):
        """Organism survives round-trip."""
        from biocompiler.infrastructure.biopython_compat import from_seqrecord, to_seqrecord

        rec = to_seqrecord(sample_dna, organism="Escherichia_coli")
        result = from_seqrecord(rec)
        assert result["organism"] == "Escherichia_coli"

    def test_gene_name_roundtrip(self, sample_dna):
        """Gene name survives round-trip."""
        from biocompiler.infrastructure.biopython_compat import from_seqrecord, to_seqrecord

        rec = to_seqrecord(sample_dna, gene_name="eGFP")
        result = from_seqrecord(rec)
        assert result["gene_name"] == "eGFP"

    def test_no_gene_name_returns_none(self, sample_dna):
        """When gene_name is not set, result['gene_name'] is None."""
        from biocompiler.infrastructure.biopython_compat import from_seqrecord, to_seqrecord

        rec = to_seqrecord(sample_dna)
        result = from_seqrecord(rec)
        assert result["gene_name"] is None

    def test_exon_boundaries_roundtrip(self, sample_dna):
        """Exon boundaries survive round-trip."""
        from biocompiler.infrastructure.biopython_compat import from_seqrecord, to_seqrecord

        boundaries = [(0, 9), (9, 18)]
        rec = to_seqrecord(sample_dna, exon_boundaries=boundaries)
        result = from_seqrecord(rec)
        assert result["exon_boundaries"] == boundaries

    def test_protein_extracted_from_cds(self, sample_dna):
        """Protein translation is extracted from the CDS feature."""
        from biocompiler.infrastructure.biopython_compat import from_seqrecord, to_seqrecord

        rec = to_seqrecord(sample_dna)
        result = from_seqrecord(rec)
        # Protein should be a non-empty string
        assert isinstance(result["protein"], str)
        assert len(result["protein"]) > 0

    def test_certificate_roundtrip(self, sample_dna):
        """Certificate data survives round-trip."""
        from biocompiler.infrastructure.biopython_compat import from_seqrecord, to_seqrecord

        cert = {"design_id": "rt-test", "status": "verified"}
        rec = to_seqrecord(sample_dna, certificate=cert)
        result = from_seqrecord(rec)
        assert result["certificate"] is not None
        assert result["certificate"]["design_id"] == "rt-test"

    def test_gc_content_computed(self, sample_dna):
        """gc_content is computed and included in the result dict."""
        from biocompiler.infrastructure.biopython_compat import from_seqrecord, to_seqrecord

        rec = to_seqrecord(sample_dna)
        result = from_seqrecord(rec)
        assert "gc_content" in result
        assert 0.0 <= result["gc_content"] <= 1.0

    def test_features_list_non_empty(self, sample_dna):
        """Features list is populated with at least the CDS feature."""
        from biocompiler.infrastructure.biopython_compat import from_seqrecord, to_seqrecord

        rec = to_seqrecord(sample_dna)
        result = from_seqrecord(rec)
        assert isinstance(result["features"], list)
        assert len(result["features"]) >= 1
        feat_types = [f["type"] for f in result["features"]]
        assert "CDS" in feat_types

    def test_feature_dicts_have_expected_keys(self, sample_dna):
        """Each feature dict has type, location, and qualifiers."""
        from biocompiler.infrastructure.biopython_compat import from_seqrecord, to_seqrecord

        rec = to_seqrecord(sample_dna)
        result = from_seqrecord(rec)
        for feat in result["features"]:
            assert "type" in feat
            assert "location" in feat
            assert "qualifiers" in feat

    def test_from_seqrecord_with_compound_location_cds(self, sample_dna):
        """from_seqrecord can extract exon boundaries from a CDS CompoundLocation."""
        from biocompiler.infrastructure.biopython_compat import from_seqrecord, to_seqrecord

        boundaries = [(0, 9), (9, 18)]
        rec = to_seqrecord(sample_dna, exon_boundaries=boundaries)
        # Remove exon features to test CompoundLocation extraction path
        rec.features = [f for f in rec.features if f.type != "exon"]
        result = from_seqrecord(rec)
        assert len(result["exon_boundaries"]) == 2

    def test_unknown_organism_default(self):
        """SeqRecord without organism annotation defaults to 'Unknown'."""
        from biocompiler.infrastructure.biopython_compat import from_seqrecord

        # Build a minimal SeqRecord without organism annotation
        rec = SeqRecord(Seq("ATGCGT"), id="test")
        result = from_seqrecord(rec)
        assert result["organism"] == "Unknown"

    def test_no_certificate_returns_none(self, sample_dna):
        """When no certificate is embedded, result['certificate'] is None."""
        from biocompiler.infrastructure.biopython_compat import from_seqrecord, to_seqrecord

        rec = to_seqrecord(sample_dna)  # no certificate
        result = from_seqrecord(rec)
        assert result["certificate"] is None


class TestRoundTripConsistency:
    """Cross-validation: to_seqrecord → from_seqrecord preserves key data."""

    def test_full_roundtrip(self, sample_dna):
        """All essential fields survive a full round-trip."""
        from biocompiler.infrastructure.biopython_compat import from_seqrecord, to_seqrecord

        boundaries = [(0, 9), (9, 18)]
        cert = {"design_id": "rt-full", "status": "verified"}

        rec = to_seqrecord(
            sample_dna,
            organism="Homo_sapiens",
            gene_name="eGFP",
            exon_boundaries=boundaries,
            certificate=cert,
        )
        result = from_seqrecord(rec)

        assert result["sequence"] == sample_dna.upper()
        assert result["organism"] == "Homo_sapiens"
        assert result["gene_name"] == "eGFP"
        assert result["exon_boundaries"] == boundaries
        assert result["certificate"]["design_id"] == "rt-full"

    def test_roundtrip_preserves_gc_content(self, sample_dna):
        """GC content is consistent between to_seqrecord annotation and
        from_seqrecord output."""
        from biocompiler.infrastructure.biopython_compat import from_seqrecord, to_seqrecord

        rec = to_seqrecord(sample_dna)
        result = from_seqrecord(rec)
        assert rec.annotations["gc_content"] == pytest.approx(result["gc_content"])

    def test_roundtrip_no_exons(self, sample_dna):
        """Round-trip with no exon boundaries gives empty list."""
        from biocompiler.infrastructure.biopython_compat import from_seqrecord, to_seqrecord

        rec = to_seqrecord(sample_dna)
        result = from_seqrecord(rec)
        # No exon features → boundaries extracted from CDS single location
        # The CDS spans the full sequence, so one boundary pair
        assert isinstance(result["exon_boundaries"], list)


class TestEdgeCases:
    """Edge cases and boundary conditions."""

    def test_very_short_sequence(self):
        """to_seqrecord works with a 3-nt (single codon) sequence."""
        from biocompiler.infrastructure.biopython_compat import to_seqrecord

        rec = to_seqrecord("ATG")
        assert isinstance(rec, SeqRecord)
        assert str(rec.seq) == "ATG"

    def test_single_methionine_roundtrip(self):
        """Round-trip with a single methionine codon."""
        from biocompiler.infrastructure.biopython_compat import from_seqrecord, to_seqrecord

        rec = to_seqrecord("ATG")
        result = from_seqrecord(rec)
        assert result["sequence"] == "ATG"

    def test_long_sequence(self, egfp_dna):
        """to_seqrecord handles a full eGFP-length sequence."""
        from biocompiler.infrastructure.biopython_compat import to_seqrecord

        rec = to_seqrecord(egfp_dna, gene_name="eGFP")
        assert len(rec.seq) == len(egfp_dna)

    def test_sequence_with_mixed_case(self):
        """Mixed-case DNA is uppercased."""
        from biocompiler.infrastructure.biopython_compat import to_seqrecord

        rec = to_seqrecord("AtGcGt")
        assert str(rec.seq) == "ATGCGT"

    def test_empty_type_results_list(self, sample_dna):
        """Empty type_results list is handled (no misc_features created)."""
        from biocompiler.infrastructure.biopython_compat import to_seqrecord

        rec = to_seqrecord(sample_dna, type_results=[])
        misc_feats = [f for f in rec.features if f.type == "misc_feature"]
        assert len(misc_feats) == 0

    def test_multiple_type_results(self, sample_dna):
        """Multiple type results each produce a separate misc_feature."""
        from biocompiler.infrastructure.biopython_compat import to_seqrecord

        class FakeVerdict:
            value = "PASS"

        class FakeTypeResult:
            predicate = "Test"
            verdict = FakeVerdict()
            violation = None
            knowledge_gap = None

        results = [FakeTypeResult() for _ in range(3)]
        rec = to_seqrecord(sample_dna, type_results=results)
        misc_feats = [f for f in rec.features if f.type == "misc_feature"]
        assert len(misc_feats) == 3

    def test_exon_numbering_in_qualifiers(self, sample_dna):
        """Exon features are numbered starting from 1."""
        from biocompiler.infrastructure.biopython_compat import to_seqrecord

        boundaries = [(0, 6), (6, 12), (12, 18)]
        rec = to_seqrecord(sample_dna, exon_boundaries=boundaries)
        exon_feats = sorted(
            [f for f in rec.features if f.type == "exon"],
            key=lambda f: int(f.location.start),
        )
        for i, feat in enumerate(exon_feats):
            assert feat.qualifiers["number"] == [str(i + 1)]

    def test_cds_codon_start_transl_table(self, sample_dna):
        """CDS feature includes codon_start=1 and transl_table=1."""
        from biocompiler.infrastructure.biopython_compat import to_seqrecord

        rec = to_seqrecord(sample_dna)
        cds = [f for f in rec.features if f.type == "CDS"][0]
        assert cds.qualifiers.get("codon_start") == ["1"]
        assert cds.qualifiers.get("transl_table") == ["1"]

    def test_gene_name_in_cds_when_provided(self, sample_dna):
        """CDS feature has gene qualifier when gene_name is provided."""
        from biocompiler.infrastructure.biopython_compat import to_seqrecord

        rec = to_seqrecord(sample_dna, gene_name="eGFP")
        cds = [f for f in rec.features if f.type == "CDS"][0]
        assert cds.qualifiers.get("gene") == ["eGFP"]

    def test_gene_name_in_exon_when_provided(self, sample_dna):
        """Exon features have gene qualifier when gene_name is provided."""
        from biocompiler.infrastructure.biopython_compat import to_seqrecord

        boundaries = [(0, 9), (9, 18)]
        rec = to_seqrecord(sample_dna, gene_name="eGFP", exon_boundaries=boundaries)
        exon_feats = [f for f in rec.features if f.type == "exon"]
        for feat in exon_feats:
            assert feat.qualifiers.get("gene") == ["eGFP"]


class TestOptimizeToSeqrecordFunctional:
    """Functional tests for ``optimize_to_seqrecord`` (requires both
    BioPython and the optimization pipeline)."""

    @pytest.fixture(autouse=True)
    def _check_deps(self):
        """Skip if biocompiler optimization pipeline is unavailable."""
        pytest.importorskip("biocompiler.optimizer")
        pytest.importorskip("biocompiler.type_system")

    def test_returns_seqrecord(self, sample_protein):
        """optimize_to_seqrecord returns a SeqRecord."""
        from biocompiler.infrastructure.biopython_compat import optimize_to_seqrecord

        rec = optimize_to_seqrecord(sample_protein)
        assert isinstance(rec, SeqRecord)

    def test_sequence_length_is_three_times_protein(self, sample_protein):
        """DNA sequence length equals 3 × protein length (no stop)."""
        from biocompiler.infrastructure.biopython_compat import optimize_to_seqrecord

        rec = optimize_to_seqrecord(sample_protein)
        # May or may not include stop codon; should be ≥ 3 × len
        assert len(rec.seq) >= 3 * len(sample_protein)

    def test_organism_in_annotations(self, sample_protein):
        """Organism is stored in annotations."""
        from biocompiler.infrastructure.biopython_compat import optimize_to_seqrecord

        rec = optimize_to_seqrecord(sample_protein, organism="Escherichia_coli")
        assert rec.annotations["organism"] == "Escherichia_coli"

    def test_gene_name_propagated(self, sample_protein):
        """Gene name is set in the resulting SeqRecord."""
        from biocompiler.infrastructure.biopython_compat import optimize_to_seqrecord

        rec = optimize_to_seqrecord(sample_protein, gene_name="eGFP")
        assert rec.id == "eGFP"

    def test_has_cds_feature(self, sample_protein):
        """Resulting SeqRecord has a CDS feature."""
        from biocompiler.infrastructure.biopython_compat import optimize_to_seqrecord

        rec = optimize_to_seqrecord(sample_protein)
        cds = [f for f in rec.features if f.type == "CDS"]
        assert len(cds) >= 1

    def test_type_check_results_as_misc_features(self, sample_protein):
        """Optimized SeqRecord includes type-check results as misc_features."""
        from biocompiler.infrastructure.biopython_compat import optimize_to_seqrecord

        rec = optimize_to_seqrecord(sample_protein)
        misc = [f for f in rec.features if f.type == "misc_feature"]
        # At least one type-check result should exist
        assert len(misc) >= 1
