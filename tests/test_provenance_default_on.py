"""
Tests for default-on provenance tracking and ProvenanceStore persistence.

Agent 3: Change provenance tracking from opt-in to default-on, and make it
persistent via ProvenanceStore.

These tests verify:
1. optimize() tracks provenance by default (track_provenance=True)
2. ProvenanceStore saves/loads correctly
3. ProvenanceStore query works
4. ProvenanceStore audit trail export works
"""

import json
import inspect
import tempfile
from pathlib import Path

import pytest

from biocompiler.optimization import optimize_sequence
from biocompiler.decision_provenance import (
    CodonDecision,
    ConstraintDecision,
    OptimizationDecisionTrail,
    DecisionProvenanceCollector,
    ProvenanceStore,
)


# ---------------------------------------------------------------------------
# Fixtures
# ---------------------------------------------------------------------------

@pytest.fixture
def store_dir(tmp_path):
    """Provide a temporary directory for ProvenanceStore tests."""
    return str(tmp_path / "provenance_store")


@pytest.fixture
def store(store_dir):
    """Provide a fresh ProvenanceStore backed by a temp directory."""
    return ProvenanceStore(store_dir=store_dir)


@pytest.fixture
def sample_trail():
    """Create a minimal OptimizationDecisionTrail for testing."""
    collector = DecisionProvenanceCollector()
    collector.start_optimization(
        protein="MVLSPADKTN",
        organism="Homo_sapiens",
        constraints=["GCInRange", "NoRestrictionSite"],
        gene_name="test_gene",
        solver_backend="greedy",
        seed=42,
    )
    collector.record_codon_decision(CodonDecision(
        position=0,
        amino_acid="M",
        original_codon=None,
        chosen_codon="ATG",
        alternatives_considered=[],
        constraint_reason="maximize_cai",
        confidence=1.0,
    ))
    collector.record_constraint_decision(ConstraintDecision(
        constraint_name="GCInRange",
        constraint_type="hard",
        action_taken="satisfied",
        positions_affected=[0, 1, 2],
        tradeoff_description="GC content within target range",
        impact_on_cai=0.0,
    ))
    return collector.finalize(output_dna="ATGGTGCTG", cai=0.95, gc=0.55)


# ---------------------------------------------------------------------------
# 1. Verify optimize() tracks provenance by default
# ---------------------------------------------------------------------------

class TestProvenanceDefaultOn:
    """Test that provenance tracking is enabled by default."""

    def test_optimize_signature_default_is_true(self):
        """The optimize_sequence function signature should have track_provenance=True."""
        sig = inspect.signature(optimize_sequence)
        param = sig.parameters["track_provenance"]
        assert param.default is True, (
            f"track_provenance default should be True, got {param.default}"
        )

    def test_optimize_produces_decision_trail_by_default(self):
        """optimize_sequence() should produce a decision_trail when called with
        default parameters (no explicit track_provenance)."""
        result = optimize_sequence("MVLSPADKTN", organism="ecoli")
        assert result.decision_trail is not None, (
            "decision_trail should not be None when track_provenance defaults to True"
        )
        assert isinstance(result.decision_trail, OptimizationDecisionTrail)

    def test_optimize_decision_trail_has_codon_decisions(self):
        """The decision trail should contain codon decisions for each amino acid."""
        result = optimize_sequence("MVLSPADKTN", organism="ecoli")
        trail = result.decision_trail
        assert trail is not None
        # 10 amino acids in "MVLSPADKTN"
        assert len(trail.codon_decisions) == 10, (
            f"Expected 10 codon decisions, got {len(trail.codon_decisions)}"
        )

    def test_optimize_explicit_false_disables_provenance(self):
        """Explicitly passing track_provenance=False should disable provenance."""
        result = optimize_sequence("MVLSPADKTN", organism="ecoli", track_provenance=False)
        assert result.decision_trail is None, (
            "decision_trail should be None when track_provenance=False"
        )

    def test_decision_trail_has_organism(self):
        """The decision trail should record the target organism."""
        result = optimize_sequence("MVLSPADKTN", organism="ecoli")
        trail = result.decision_trail
        assert trail is not None
        assert trail.organism == "Escherichia_coli"


# ---------------------------------------------------------------------------
# 2. ProvenanceStore save/load
# ---------------------------------------------------------------------------

class TestProvenanceStoreSaveLoad:
    """Test ProvenanceStore save and load functionality."""

    def test_save_returns_record_id(self, store, sample_trail):
        """save() should return a string record ID."""
        record_id = store.save(sample_trail)
        assert isinstance(record_id, str)
        assert len(record_id) > 0

    def test_save_creates_json_file(self, store, sample_trail):
        """save() should create a JSON file in the store directory."""
        record_id = store.save(sample_trail)
        filepath = store.store_dir / f"{record_id}.json"
        assert filepath.exists()

    def test_load_roundtrip(self, store, sample_trail):
        """Loading a saved record should produce an equivalent trail."""
        record_id = store.save(sample_trail)
        loaded = store.load(record_id)
        assert loaded.gene_name == sample_trail.gene_name
        assert loaded.input_protein == sample_trail.input_protein
        assert loaded.output_dna == sample_trail.output_dna
        assert loaded.organism == sample_trail.organism
        assert loaded.total_cai == sample_trail.total_cai
        assert loaded.total_gc == sample_trail.total_gc
        assert len(loaded.codon_decisions) == len(sample_trail.codon_decisions)
        assert len(loaded.constraint_decisions) == len(sample_trail.constraint_decisions)

    def test_load_nonexistent_raises_file_not_found(self, store):
        """load() with an invalid ID should raise FileNotFoundError."""
        with pytest.raises(FileNotFoundError, match="not found"):
            store.load("nonexistent-record-id")

    def test_load_corrupted_file_raises_value_error(self, store):
        """load() with a corrupted JSON file should raise ValueError."""
        filepath = store.store_dir / "corrupted.json"
        filepath.write_text("NOT VALID JSON {{{", encoding="utf-8")
        with pytest.raises(ValueError, match="Corrupted"):
            store.load("corrupted")

    def test_save_multiple_records(self, store, sample_trail):
        """Saving multiple records should produce different IDs."""
        id1 = store.save(sample_trail)
        id2 = store.save(sample_trail)
        assert id1 != id2

    def test_store_dir_created_automatically(self, tmp_path):
        """ProvenanceStore should create the store directory if it doesn't exist."""
        new_dir = str(tmp_path / "auto_created" / "provenance")
        assert not Path(new_dir).exists()
        store = ProvenanceStore(store_dir=new_dir)
        assert Path(new_dir).exists()
        assert store.store_dir == Path(new_dir)


# ---------------------------------------------------------------------------
# 3. ProvenanceStore query
# ---------------------------------------------------------------------------

class TestProvenanceStoreQuery:
    """Test ProvenanceStore query functionality."""

    def test_query_by_organism(self, store, sample_trail):
        """Querying by organism should return only matching records."""
        store.save(sample_trail)
        # Create another trail with different organism
        collector2 = DecisionProvenanceCollector()
        collector2.start_optimization(
            protein="MVLSPADKTN",
            organism="Saccharomyces_cerevisiae",
            constraints=[],
            gene_name="test_gene",
            solver_backend="greedy",
        )
        trail2 = collector2.finalize(output_dna="ATGGTGCTG", cai=0.80, gc=0.40)
        store.save(trail2)

        results = store.query(organism="Homo_sapiens")
        assert len(results) == 1
        assert results[0].organism == "Homo_sapiens"

    def test_query_by_protein_name(self, store, sample_trail):
        """Querying by protein_name (gene_name) should return only matching records."""
        store.save(sample_trail)
        # Create a trail without gene_name
        collector2 = DecisionProvenanceCollector()
        collector2.start_optimization(
            protein="MVLSPADKTN",
            organism="Homo_sapiens",
            constraints=[],
            solver_backend="greedy",
        )
        trail2 = collector2.finalize(output_dna="ATGGTGCTG", cai=0.80, gc=0.40)
        store.save(trail2)

        results = store.query(protein_name="test_gene")
        assert len(results) == 1
        assert results[0].gene_name == "test_gene"

    def test_query_no_filters_returns_all(self, store, sample_trail):
        """Query with no filters should return all stored records."""
        store.save(sample_trail)
        store.save(sample_trail)
        results = store.query()
        assert len(results) == 2

    def test_query_by_date_range(self, store, sample_trail):
        """Querying by date_range should filter by timestamp."""
        store.save(sample_trail)
        # Wide date range should match
        results = store.query(
            date_range=("2000-01-01T00:00:00+00:00", "2100-01-01T00:00:00+00:00")
        )
        assert len(results) >= 1

        # Narrow past range should not match
        results = store.query(
            date_range=("2000-01-01T00:00:00+00:00", "2001-01-01T00:00:00+00:00")
        )
        assert len(results) == 0

    def test_query_combined_filters(self, store, sample_trail):
        """Multiple filters should be combined with AND logic."""
        store.save(sample_trail)
        # This should match: correct organism AND wide date range
        results = store.query(
            organism="Homo_sapiens",
            date_range=("2000-01-01T00:00:00+00:00", "2100-01-01T00:00:00+00:00"),
        )
        assert len(results) == 1

        # This should not match: wrong organism
        results = store.query(
            organism="Escherichia_coli",
            date_range=("2000-01-01T00:00:00+00:00", "2100-01-01T00:00:00+00:00"),
        )
        assert len(results) == 0


# ---------------------------------------------------------------------------
# 4. ProvenanceStore audit trail export
# ---------------------------------------------------------------------------

class TestProvenanceStoreAuditTrail:
    """Test ProvenanceStore export_audit_trail functionality."""

    def test_export_json_format(self, store, sample_trail):
        """Exporting in JSON format should produce valid JSON."""
        record_id = store.save(sample_trail)
        audit_json = store.export_audit_trail(record_id, format="json")
        parsed = json.loads(audit_json)
        assert parsed["organism"] == "Homo_sapiens"
        assert parsed["gene_name"] == "test_gene"
        assert len(parsed["codon_decisions"]) == 1

    def test_export_text_format(self, store, sample_trail):
        """Exporting in text format should produce a human-readable summary."""
        record_id = store.save(sample_trail)
        audit_text = store.export_audit_trail(record_id, format="text")
        assert "Homo_sapiens" in audit_text
        assert "test_gene" in audit_text
        assert "Codon decisions" in audit_text
        assert "M -> ATG" in audit_text

    def test_export_unsupported_format_raises(self, store, sample_trail):
        """Exporting in an unsupported format should raise ValueError."""
        record_id = store.save(sample_trail)
        with pytest.raises(ValueError, match="Unsupported audit trail format"):
            store.export_audit_trail(record_id, format="xml")

    def test_export_nonexistent_record_raises(self, store):
        """Exporting a nonexistent record should raise FileNotFoundError."""
        with pytest.raises(FileNotFoundError):
            store.export_audit_trail("nonexistent-id")


# ---------------------------------------------------------------------------
# 5. Integration: optimize_sequence result through ProvenanceStore
# ---------------------------------------------------------------------------

class TestProvenanceIntegration:
    """Integration test: full pipeline from optimization through persistent storage."""

    def test_optimize_save_load_roundtrip(self, store):
        """Optimize a protein, save the trail, load it back, verify consistency."""
        result = optimize_sequence("MVLSPADKTN", organism="ecoli")
        trail = result.decision_trail
        assert trail is not None

        record_id = store.save(trail)
        loaded = store.load(record_id)

        assert loaded.input_protein == trail.input_protein
        assert loaded.output_dna == trail.output_dna
        assert loaded.organism == trail.organism
        assert loaded.total_cai == trail.total_cai
        assert loaded.total_gc == trail.total_gc
        assert len(loaded.codon_decisions) == len(trail.codon_decisions)

    def test_optimize_query_by_organism(self, store):
        """Optimize for different organisms, then query by organism."""
        result_ecoli = optimize_sequence("MVLSPADKTN", organism="ecoli")
        result_human = optimize_sequence("MVLSPADKTN", organism="human")

        store.save(result_ecoli.decision_trail)
        store.save(result_human.decision_trail)

        ecoli_results = store.query(organism="Escherichia_coli")
        human_results = store.query(organism="Homo_sapiens")

        assert len(ecoli_results) == 1
        assert len(human_results) == 1
        assert ecoli_results[0].organism == "Escherichia_coli"
        assert human_results[0].organism == "Homo_sapiens"
