"""
Tests for output provenance: Certificate dataclass + export integration.

Covers:
1. Certificate dataclass construction and field access
2. compute_certificate returns CertLevel (GOLD/SILVER/BRONZE)
3. export_genbank_with_certificate includes certificate data
4. FASTA export includes GC content and organism
5. GenBank export has proper LOCUS/FEATURES/ORIGIN structure
6. Decision audit: optimization results include metadata
   (DecisionRecord, ProvenanceTracker, OptimizationProvenance)
"""

import hashlib
import json
from datetime import datetime, timezone

import pytest

from biocompiler.certificate import compute_certificate, generate_certificate
from biocompiler.export import (
    export_fasta,
    export_genbank,
    export_genbank_with_certificate,
)
from biocompiler.provenance import (
    DecisionRecord,
    OptimizationProvenance,
    ProvenanceTracker,
)
from biocompiler.type_system import CertLevel, PredicateResult, SpliceVerdict
from biocompiler.types import Certificate, TypeCheckResult, Verdict, SLOTMode


# ─── Shared fixtures ────────────────────────────────────────────────────────

SAMPLE_SEQ = "ATGGTGAGCAAGGGCGAGGAGCTGTTCACCGGGGTGGTGCCC"  # 45 bp


def _make_type_results(
    verdicts: list[tuple[str, Verdict]] | None = None,
) -> list[TypeCheckResult]:
    """Build a list of TypeCheckResult objects for testing."""
    if verdicts is None:
        verdicts = [
            ("GCInRange", Verdict.PASS),
            ("NoStopCodons", Verdict.PASS),
        ]
    return [
        TypeCheckResult(predicate=p, verdict=v, violation=None if v == Verdict.PASS else "See details")
        for p, v in verdicts
    ]


def _make_predicate_results(
    all_pass: bool = True,
    with_mutagenesis: bool = False,
    with_unavoidable: bool = False,
) -> list[PredicateResult]:
    """Build a list of PredicateResult objects for certificate level tests."""
    results = [
        PredicateResult("NoStopCodons", True, verdict=Verdict.PASS, details="No internal stop codons"),
        PredicateResult(
            "NoCrypticSplice", True, verdict=SpliceVerdict.PASS,
            details="No GT dinucleotides found",
        ),
        PredicateResult("NoCpGIsland", True, details="Worst CpG Obs/Exp ratio 0.300 <= 0.6"),
        PredicateResult("NoRestrictionSite", True, details="No restriction sites found"),
    ]
    if with_unavoidable:
        results.append(PredicateResult(
            "NoGTDinucleotide", True,
            details="All 2 GT dinucleotides are unavoidable",
        ))
    elif with_mutagenesis:
        results.append(PredicateResult(
            "NoGTDinucleotide", True,
            details="No GT dinucleotides found mutagenesis applied: pos 3:V→I",
        ))
    else:
        results.append(PredicateResult(
            "NoGTDinucleotide", True, details="No GT dinucleotides found",
        ))

    if not all_pass:
        results.append(PredicateResult(
            "ValidCodingSeq", False,
            details="Sequence length not divisible by 3",
        ))
    else:
        results.append(PredicateResult("ValidCodingSeq", True, details="All codons valid"))

    results.append(PredicateResult(
        "ConservationScore", True, details="All AA conservation scores >= -1",
    ))
    results.append(PredicateResult(
        "CodonOptimality", True, details="Worst CAI: GCT=0.7244, min=0.0",
    ))
    return results


def _make_certificate(
    design_id: str = "TEST_CERT_001",
    types: list[dict] | None = None,
    provenance: dict | None = None,
) -> Certificate:
    """Helper to construct a Certificate for export tests."""
    if types is None:
        types = [
            {"predicate": "GCInRange", "verdict": "PASS"},
            {"predicate": "NoStopCodons", "verdict": "PASS"},
        ]
    if provenance is None:
        provenance = {
            "tool": "BioCompiler",
            "version": "7.2.0",
            "timestamp": "2025-06-15T12:00:00+00:00",
            "input_hash": "abc123",
            "overall_status": "FULL_PASS",
        }
    return Certificate(
        version="7.2.0",
        design_id=design_id,
        sequence=SAMPLE_SEQ,
        types=types,
        provenance=provenance,
    )


# ═══════════════════════════════════════════════════════════════════════════════
# 1. Certificate dataclass construction and field access
# ═══════════════════════════════════════════════════════════════════════════════

class TestCertificateDataclass:
    """Tests for the Certificate dataclass: construction, fields, serialization."""

    def test_construction_with_all_fields(self):
        """Certificate can be constructed with all required fields."""
        cert = Certificate(
            version="1.0",
            design_id="ABC123",
            sequence="ATGC",
            types=[{"predicate": "gc", "verdict": "PASS"}],
            provenance={"tool": "BioCompiler", "timestamp": "2025-01-01"},
        )
        assert cert.version == "1.0"
        assert cert.design_id == "ABC123"
        assert cert.sequence == "ATGC"
        assert len(cert.types) == 1
        assert cert.provenance["tool"] == "BioCompiler"

    def test_field_access_version(self):
        cert = _make_certificate()
        assert isinstance(cert.version, str)
        assert cert.version == "7.2.0"

    def test_field_access_design_id(self):
        cert = _make_certificate(design_id="MY_DESIGN")
        assert cert.design_id == "MY_DESIGN"

    def test_field_access_sequence(self):
        cert = _make_certificate()
        assert cert.sequence == SAMPLE_SEQ

    def test_field_access_types(self):
        types = [
            {"predicate": "GCInRange", "verdict": "PASS"},
            {"predicate": "NoStopCodons", "verdict": "FAIL"},
        ]
        cert = _make_certificate(types=types)
        assert len(cert.types) == 2
        assert cert.types[1]["verdict"] == "FAIL"

    def test_field_access_provenance(self):
        prov = {"tool": "BioCompiler", "timestamp": "2025-03-01", "input_hash": "deadbeef"}
        cert = _make_certificate(provenance=prov)
        assert cert.provenance["input_hash"] == "deadbeef"
        assert cert.provenance["tool"] == "BioCompiler"

    def test_to_dict_roundtrip(self):
        """to_dict should produce a JSON-compatible dict with all fields."""
        cert = _make_certificate()
        d = cert.to_dict()
        # v2 certificates include hash_version
        expected_keys = {"version", "design_id", "sequence", "types", "provenance", "hash_version"}
        assert set(d.keys()) == expected_keys
        assert d["design_id"] == "TEST_CERT_001"
        assert d["hash_version"] == 2

    def test_from_dict_validates_required_keys(self):
        """from_dict raises ValueError if required keys are missing."""
        with pytest.raises(ValueError, match="missing keys"):
            Certificate.from_dict({"version": "1.0"})  # missing design_id, sequence, etc.

    def test_from_dict_roundtrip(self):
        """from_dict(to_dict()) should produce an equivalent Certificate."""
        cert = _make_certificate()
        d = cert.to_dict()
        restored = Certificate.from_dict(d)
        assert restored.version == cert.version
        assert restored.design_id == cert.design_id
        assert restored.sequence == cert.sequence
        assert restored.types == cert.types
        assert restored.provenance == cert.provenance

    def test_empty_types_list(self):
        """Certificate with empty types list is valid."""
        cert = _make_certificate(types=[])
        assert cert.types == []

    def test_provenance_nested_structure(self):
        """Provenance can contain nested dicts (e.g., parameters, mutagenesis)."""
        prov = {
            "tool": "BioCompiler",
            "timestamp": "2025-01-01",
            "input_hash": "abc",
            "parameters": {"gc_lo": 0.3, "gc_hi": 0.7},
            "mutagenesis": {"applied": False},
        }
        cert = _make_certificate(provenance=prov)
        assert cert.provenance["parameters"]["gc_lo"] == 0.3
        assert cert.provenance["mutagenesis"]["applied"] is False

    def test_to_dict_is_json_serializable(self):
        """to_dict output should be serializable to JSON without errors."""
        cert = _make_certificate()
        d = cert.to_dict()
        json_str = json.dumps(d)
        assert isinstance(json_str, str)
        restored = json.loads(json_str)
        assert restored["design_id"] == "TEST_CERT_001"


# ═══════════════════════════════════════════════════════════════════════════════
# 2. compute_certificate returns CertLevel (GOLD / SILVER / BRONZE)
# ═══════════════════════════════════════════════════════════════════════════════

class TestComputeCertificate:
    """Tests for certificate level computation from predicate results."""

    def test_returns_certlevel_type(self):
        """compute_certificate should return a CertLevel enum member."""
        results = _make_predicate_results(all_pass=True)
        level = compute_certificate(results)
        assert isinstance(level, CertLevel)

    def test_gold_when_all_pass_no_mutagenesis(self):
        """All predicates pass with no mutagenesis/unavoidable → GOLD."""
        results = _make_predicate_results(all_pass=True)
        assert compute_certificate(results) == CertLevel.GOLD

    def test_silver_when_mutagenesis(self):
        """All pass but mutagenesis mentioned in details → SILVER."""
        results = _make_predicate_results(all_pass=True, with_mutagenesis=True)
        assert compute_certificate(results) == CertLevel.SILVER

    def test_silver_when_unavoidable(self):
        """All pass but unavoidable mentioned in details → SILVER."""
        results = _make_predicate_results(all_pass=True, with_unavoidable=True)
        assert compute_certificate(results) == CertLevel.SILVER

    def test_bronze_when_unsatisfied(self):
        """Any unsatisfied predicate → BRONZE."""
        results = _make_predicate_results(all_pass=False)
        assert compute_certificate(results) == CertLevel.BRONZE

    def test_generate_certificate_returns_certificate(self):
        """generate_certificate should return a Certificate dataclass."""
        type_results = _make_type_results()
        cert = generate_certificate(
            sequence=SAMPLE_SEQ,
            type_results=type_results,
            input_params={"gene": "eGFP"},
        )
        assert isinstance(cert, Certificate)

    def test_generate_certificate_has_provenance(self):
        """Generated certificate must have a provenance dict with required keys."""
        type_results = _make_type_results()
        cert = generate_certificate(
            sequence=SAMPLE_SEQ,
            type_results=type_results,
            input_params={"gene": "eGFP"},
        )
        assert "tool" in cert.provenance
        assert "version" in cert.provenance
        assert "timestamp" in cert.provenance
        assert "input_hash" in cert.provenance

    def test_generate_certificate_design_id_is_v2_hash(self):
        """design_id should be the v2 hash covering sequence + predicates + params."""
        from biocompiler.certificate import _compute_certificate_hash
        type_results = _make_type_results()
        cert = generate_certificate(
            sequence=SAMPLE_SEQ,
            type_results=type_results,
            input_params={"gene": "eGFP"},
        )
        # v1 hash (sequence-only) should NOT match
        v1_hash = hashlib.sha256(SAMPLE_SEQ.encode()).hexdigest()
        assert cert.design_id != v1_hash, "v2 hash must differ from v1 sequence-only hash"
        # v2 hash should match
        params = cert.provenance.get("parameters", {})
        expected = _compute_certificate_hash(
            sequence=SAMPLE_SEQ,
            types_list=cert.types,
            params=params,
            hash_version=2,
        )
        assert cert.design_id == expected

    def test_generate_certificate_empty_sequence_raises(self):
        """Empty sequence should raise ValueError."""
        with pytest.raises(ValueError, match="Sequence must not be empty"):
            generate_certificate("", _make_type_results(), input_params={})

    def test_generate_certificate_empty_type_results_raises(self):
        """Empty type_results should raise ValueError."""
        with pytest.raises(ValueError, match="Type results must not be empty"):
            generate_certificate(SAMPLE_SEQ, [], input_params={})

    def test_generate_certificate_graduated_mode(self):
        """In graduated mode (default), certificate is generated even with failures."""
        type_results = _make_type_results([
            ("GCInRange", Verdict.PASS),
            ("NoStopCodons", Verdict.FAIL),
        ])
        cert = generate_certificate(
            sequence=SAMPLE_SEQ,
            type_results=type_results,
            input_params={"gene": "eGFP"},
            require_all_pass=False,
        )
        assert isinstance(cert, Certificate)
        assert "PARTIAL" in cert.provenance["overall_status"]


# ═══════════════════════════════════════════════════════════════════════════════
# 3. export_genbank_with_certificate includes certificate data
# ═══════════════════════════════════════════════════════════════════════════════

class TestExportGenbankWithCertificate:
    """Tests that GenBank export with certificate includes certificate provenance."""

    def test_includes_certificate_id_in_comment(self):
        cert = _make_certificate(design_id="CERT_XYZ_789")
        result = export_genbank_with_certificate(SAMPLE_SEQ, certificate=cert)
        assert "Certificate ID:" in result

    def test_includes_certificate_timestamp_in_comment(self):
        cert = _make_certificate(provenance={
            "tool": "BioCompiler",
            "version": "7.2.0",
            "timestamp": "2025-06-15T12:00:00+00:00",
            "input_hash": "abc123",
        })
        result = export_genbank_with_certificate(SAMPLE_SEQ, certificate=cert)
        assert "Certificate timestamp:" in result
        assert "2025-06-15" in result

    def test_uses_certificate_design_id_as_locus(self):
        cert = _make_certificate(design_id="MYDESIGN123456")
        result = export_genbank_with_certificate(SAMPLE_SEQ, certificate=cert)
        locus_line = result.split("\n")[0]
        assert "MYDESIGN123456"[:16].upper() in locus_line

    def test_type_results_reconstructed_from_certificate(self):
        """Certificate types are reconstructed into TypeCheckResults for annotation."""
        types = [
            {"predicate": "GCInRange", "verdict": "PASS"},
            {"predicate": "NoStopCodons", "verdict": "FAIL"},
        ]
        cert = _make_certificate(types=types)
        result = export_genbank_with_certificate(SAMPLE_SEQ, certificate=cert)
        # FAIL verdicts should produce TYPE FAIL annotations
        assert "TYPE FAIL: NoStopCodons" in result

    def test_definition_includes_design_id_prefix(self):
        cert = _make_certificate(design_id="MYDESN_XYZ")
        result = export_genbank_with_certificate(SAMPLE_SEQ, certificate=cert)
        assert "DEFINITION" in result
        assert "MYDESN_X"[:8] in result

    def test_gene_name_propagated(self):
        cert = _make_certificate()
        result = export_genbank_with_certificate(
            SAMPLE_SEQ, certificate=cert, gene_name="eGFP",
        )
        assert '/gene="eGFP"' in result

    def test_exon_boundaries_propagated(self):
        cert = _make_certificate()
        boundaries = [(0, 20), (25, 45)]
        result = export_genbank_with_certificate(
            SAMPLE_SEQ, certificate=cert, exon_boundaries=boundaries,
        )
        assert "join(" in result
        assert "1..20" in result

    def test_valid_genbank_structure_with_certificate(self):
        """GenBank with certificate still has LOCUS/FEATURES/ORIGIN/terminator."""
        cert = _make_certificate()
        result = export_genbank_with_certificate(SAMPLE_SEQ, certificate=cert)
        assert "LOCUS" in result
        assert "FEATURES" in result
        assert "ORIGIN" in result
        assert result.rstrip().endswith("//")


# ═══════════════════════════════════════════════════════════════════════════════
# 4. FASTA export includes GC content and organism
# ═══════════════════════════════════════════════════════════════════════════════

class TestFastaExportProvenance:
    """Tests that FASTA export includes GC content and organism metadata."""

    def test_header_includes_gc_content(self):
        result = export_fasta(SAMPLE_SEQ)
        header = result.split("\n")[0]
        gc_parts = [p for p in header.split("|") if p.startswith("gc=")]
        assert len(gc_parts) == 1
        gc_val = float(gc_parts[0].split("=")[1])
        assert 0.0 <= gc_val <= 1.0

    def test_gc_content_value_correct_all_gc(self):
        """Sequence of all G/C should report gc=1.000."""
        result = export_fasta("GCGCGCGC")
        header = result.split("\n")[0]
        for part in header.split("|"):
            if part.startswith("gc="):
                assert float(part.split("=")[1]) == 1.0

    def test_gc_content_value_correct_all_at(self):
        """Sequence of all A/T should report gc=0.000."""
        result = export_fasta("ATATATAT")
        header = result.split("\n")[0]
        for part in header.split("|"):
            if part.startswith("gc="):
                assert float(part.split("=")[1]) == 0.0

    def test_gc_content_mixed(self):
        """ATGC should report gc=0.5."""
        result = export_fasta("ATGC")
        header = result.split("\n")[0]
        for part in header.split("|"):
            if part.startswith("gc="):
                assert float(part.split("=")[1]) == 0.5

    def test_header_includes_organism(self):
        result = export_fasta(SAMPLE_SEQ, organism="Homo_sapiens")
        assert "organism=Homo_sapiens" in result

    def test_organism_custom(self):
        result = export_fasta(SAMPLE_SEQ, organism="Escherichia_coli")
        assert "organism=Escherichia_coli" in result

    def test_header_includes_length(self):
        result = export_fasta(SAMPLE_SEQ)
        header = result.split("\n")[0]
        for part in header.split("|"):
            if part.startswith("len="):
                assert str(len(SAMPLE_SEQ)) in part

    def test_header_includes_protein_len(self):
        result = export_fasta(SAMPLE_SEQ)
        assert "protein_len=" in result
        assert "aa" in result

    def test_gc_format_three_decimals(self):
        """GC content should be formatted to 3 decimal places."""
        result = export_fasta("ATGCATGCATGC")
        header = result.split("\n")[0]
        for part in header.split("|"):
            if part.startswith("gc="):
                val = part.split("=")[1]
                assert "." in val
                assert len(val.split(".")[1]) == 3


# ═══════════════════════════════════════════════════════════════════════════════
# 5. GenBank export has proper LOCUS/FEATURES/ORIGIN structure
# ═══════════════════════════════════════════════════════════════════════════════

class TestGenBankStructure:
    """Tests that GenBank output has correct structural ordering."""

    def test_has_locus(self):
        result = export_genbank(SAMPLE_SEQ)
        assert "LOCUS" in result

    def test_has_features(self):
        result = export_genbank(SAMPLE_SEQ)
        assert "FEATURES" in result

    def test_has_origin(self):
        result = export_genbank(SAMPLE_SEQ)
        assert "ORIGIN" in result

    def test_has_terminator(self):
        result = export_genbank(SAMPLE_SEQ)
        assert result.rstrip().endswith("//")

    def test_locus_before_features(self):
        result = export_genbank(SAMPLE_SEQ)
        assert result.index("LOCUS") < result.index("FEATURES")

    def test_features_before_origin(self):
        result = export_genbank(SAMPLE_SEQ)
        assert result.index("FEATURES") < result.index("ORIGIN")

    def test_origin_before_terminator(self):
        result = export_genbank(SAMPLE_SEQ)
        assert result.index("ORIGIN") < result.index("//")

    def test_locus_line_format(self):
        """LOCUS line should contain name, bp count, molecule type, topology, division."""
        result = export_genbank(
            SAMPLE_SEQ, locus_name="MYGENE", molecule_type="DNA", topology="linear",
        )
        locus_line = result.split("\n")[0]
        assert "MYGENE" in locus_line
        assert "bp" in locus_line
        assert "DNA" in locus_line
        assert "linear" in locus_line
        assert "SYN" in locus_line

    def test_locus_includes_sequence_length(self):
        result = export_genbank(SAMPLE_SEQ)
        assert f"{len(SAMPLE_SEQ)} bp" in result

    def test_has_definition_after_locus(self):
        result = export_genbank(SAMPLE_SEQ)
        assert result.index("LOCUS") < result.index("DEFINITION")

    def test_has_accession(self):
        result = export_genbank(SAMPLE_SEQ)
        assert "ACCESSION" in result

    def test_accession_before_version(self):
        result = export_genbank(SAMPLE_SEQ)
        assert result.index("ACCESSION") < result.index("VERSION")

    def test_has_source_and_organism(self):
        result = export_genbank(SAMPLE_SEQ, organism="Homo_sapiens")
        assert "SOURCE" in result
        assert "ORGANISM" in result
        assert "Homo_sapiens" in result

    def test_features_table_has_gene_when_name_provided(self):
        result = export_genbank(SAMPLE_SEQ, gene_name="eGFP")
        assert '/gene="eGFP"' in result

    def test_features_has_cds(self):
        result = export_genbank(SAMPLE_SEQ)
        assert "CDS" in result

    def test_origin_contains_numbered_sequence(self):
        """ORIGIN section should contain the sequence in numbered format."""
        result = export_genbank("ATGCATGCATGC")
        origin_start = result.index("ORIGIN")
        terminator = result.index("//")
        origin_section = result[origin_start:terminator]
        # Should start with a number
        import re
        lines = origin_section.split("\n")[1:]  # skip "ORIGIN" line
        for line in lines:
            line = line.strip()
            if line:
                assert re.match(r"^\d+", line), f"Line doesn't start with number: {line}"

    def test_terminator_is_last_non_whitespace(self):
        result = export_genbank(SAMPLE_SEQ)
        assert result.rstrip().endswith("//")


# ═══════════════════════════════════════════════════════════════════════════════
# 6. Decision audit: optimization results include metadata
# ═══════════════════════════════════════════════════════════════════════════════

class TestDecisionRecord:
    """Tests for DecisionRecord dataclass construction and serialization."""

    def test_construction(self):
        rec = DecisionRecord(
            timestamp="2025-06-15T12:00:00+00:00",
            decision_type="codon_selected",
            position=12,
            chosen_value="GTC",
            alternatives_considered=["GTG", "GTA", "GTT"],
            rationale="Highest CAI codon for Valine",
            constraint_context={"cai": 0.92, "gc": 0.54},
        )
        assert rec.decision_type == "codon_selected"
        assert rec.position == 12
        assert rec.chosen_value == "GTC"
        assert len(rec.alternatives_considered) == 3

    def test_to_dict(self):
        rec = DecisionRecord(
            timestamp="2025-06-15T12:00:00+00:00",
            decision_type="codon_selected",
            position=0,
            chosen_value="ATG",
            alternatives_considered=[],
            rationale="Start codon",
            constraint_context={},
        )
        d = rec.to_dict()
        assert d["decision_type"] == "codon_selected"
        assert d["chosen_value"] == "ATG"
        assert "timestamp" in d

    def test_from_dict_roundtrip(self):
        rec = DecisionRecord(
            timestamp="2025-06-15T12:00:00+00:00",
            decision_type="mutation_applied",
            position=42,
            chosen_value="V42I",
            alternatives_considered=["V42L"],
            rationale="Remove GT dinucleotide",
            constraint_context={"blosum62": 3},
        )
        d = rec.to_dict()
        restored = DecisionRecord.from_dict(d)
        assert restored.decision_type == rec.decision_type
        assert restored.position == rec.position
        assert restored.chosen_value == rec.chosen_value
        assert restored.alternatives_considered == rec.alternatives_considered
        assert restored.rationale == rec.rationale
        assert restored.constraint_context == rec.constraint_context

    def test_from_dict_missing_keys_raises(self):
        with pytest.raises(ValueError, match="missing keys"):
            DecisionRecord.from_dict({"timestamp": "now", "decision_type": "x"})

    def test_to_dict_is_json_serializable(self):
        rec = DecisionRecord(
            timestamp="2025-06-15T12:00:00+00:00",
            decision_type="codon_selected",
            position=0,
            chosen_value="ATG",
            alternatives_considered=[],
            rationale="Start",
            constraint_context={"cai": 0.95},
        )
        json_str = json.dumps(rec.to_dict())
        assert isinstance(json_str, str)


class TestProvenanceTracker:
    """Tests for ProvenanceTracker: recording, querying, serialization."""

    def test_record_and_retrieve_decision(self):
        tracker = ProvenanceTracker(seed=42)
        rec = DecisionRecord(
            timestamp=datetime.now(timezone.utc).isoformat(),
            decision_type="codon_selected",
            position=12,
            chosen_value="GTC",
            alternatives_considered=["GTG"],
            rationale="Best CAI",
            constraint_context={"cai": 0.9},
        )
        tracker.record_decision(rec)
        decisions = tracker.get_decisions_for_position(12)
        assert len(decisions) == 1
        assert decisions[0].chosen_value == "GTC"

    def test_get_decisions_for_empty_position(self):
        tracker = ProvenanceTracker(seed=0)
        assert tracker.get_decisions_for_position(999) == []

    def test_multiple_decisions_at_same_position(self):
        tracker = ProvenanceTracker(seed=0)
        for val in ["GTC", "GTG"]:
            tracker.record_decision(DecisionRecord(
                timestamp="2025-01-01T00:00:00Z",
                decision_type="codon_selected",
                position=5,
                chosen_value=val,
                alternatives_considered=[],
                rationale="",
                constraint_context={},
            ))
        decisions = tracker.get_decisions_for_position(5)
        assert len(decisions) == 2

    def test_full_audit_trail(self):
        tracker = ProvenanceTracker(seed=0)
        for i in range(5):
            tracker.record_decision(DecisionRecord(
                timestamp="2025-01-01T00:00:00Z",
                decision_type="codon_selected",
                position=i * 3,
                chosen_value="ATG",
                alternatives_considered=[],
                rationale="",
                constraint_context={},
            ))
        trail = tracker.get_full_audit_trail()
        assert len(trail) == 5

    def test_len(self):
        tracker = ProvenanceTracker(seed=0)
        assert len(tracker) == 0
        tracker.record_decision(DecisionRecord(
            timestamp="2025-01-01", decision_type="x", position=0,
            chosen_value="A", alternatives_considered=[], rationale="",
            constraint_context={},
        ))
        assert len(tracker) == 1

    def test_seed_stored(self):
        tracker = ProvenanceTracker(seed=12345)
        assert tracker.seed == 12345

    def test_to_dict_includes_seed(self):
        tracker = ProvenanceTracker(seed=99)
        d = tracker.to_dict()
        assert d["seed"] == 99
        assert "decision_count" in d
        assert "decisions" in d

    def test_to_dict_and_from_dict_roundtrip(self):
        tracker = ProvenanceTracker(seed=42)
        tracker.record_decision(DecisionRecord(
            timestamp="2025-06-15T12:00:00+00:00",
            decision_type="codon_selected",
            position=6,
            chosen_value="GCT",
            alternatives_considered=["GCC", "GCA"],
            rationale="Good CAI",
            constraint_context={"gc": 0.55},
        ))
        d = tracker.to_dict()
        restored = ProvenanceTracker.from_dict(d)
        assert restored.seed == 42
        assert len(restored) == 1
        decisions = restored.get_decisions_for_position(6)
        assert decisions[0].chosen_value == "GCT"

    def test_to_json_and_from_json_roundtrip(self):
        tracker = ProvenanceTracker(seed=7)
        tracker.record_decision(DecisionRecord(
            timestamp="2025-01-01", decision_type="mutation_applied",
            position=10, chosen_value="V10I",
            alternatives_considered=[], rationale="BLOSUM62>=0",
            constraint_context={"blosum62": 3},
        ))
        json_str = tracker.to_json()
        restored = ProvenanceTracker.from_json(json_str)
        assert restored.seed == 7
        assert len(restored) == 1

    def test_from_dict_missing_seed_raises(self):
        with pytest.raises(ValueError, match="missing key 'seed'"):
            ProvenanceTracker.from_dict({})

    def test_record_non_decision_record_raises_type_error(self):
        tracker = ProvenanceTracker(seed=0)
        with pytest.raises(TypeError, match="Expected DecisionRecord"):
            tracker.record_decision("not a record")  # type: ignore[arg-type]

    def test_repr(self):
        tracker = ProvenanceTracker(seed=42)
        assert "seed=42" in repr(tracker)
        assert "decisions=0" in repr(tracker)


class TestOptimizationProvenance:
    """Tests for OptimizationProvenance: end-to-end provenance snapshot."""

    def _make_provenance(self) -> OptimizationProvenance:
        return OptimizationProvenance(
            input_protein="MVLSPADKTNVKAAWGKVGAHAGEYGAEALERMFLSFPTTKTYFPHFDLSH",
            organism="Homo_sapiens",
            solver_backend="greedy",
            config_snapshot={"gc_lo": 0.30, "gc_hi": 0.70, "cai_threshold": 0.5},
            decisions=[
                DecisionRecord(
                    timestamp="2025-06-15T12:00:00+00:00",
                    decision_type="codon_selected",
                    position=0,
                    chosen_value="ATG",
                    alternatives_considered=[],
                    rationale="Start codon, fixed",
                    constraint_context={},
                ),
            ],
            final_sequence="ATGGTGCTGCCTGCTCCTGCTG",
            solve_time_seconds=1.23,
            constraints_active=["GCInRange", "NoCrypticSplice", "NoRestrictionSite"],
        )

    def test_construction(self):
        prov = self._make_provenance()
        assert prov.organism == "Homo_sapiens"
        assert prov.solver_backend == "greedy"
        assert len(prov.decisions) == 1
        assert prov.solve_time_seconds == 1.23

    def test_to_dict(self):
        prov = self._make_provenance()
        d = prov.to_dict()
        assert d["organism"] == "Homo_sapiens"
        assert d["solver_backend"] == "greedy"
        assert "decisions" in d
        assert len(d["decisions"]) == 1
        assert d["solve_time_seconds"] == 1.23
        assert d["constraints_active"] == ["GCInRange", "NoCrypticSplice", "NoRestrictionSite"]

    def test_from_dict_roundtrip(self):
        prov = self._make_provenance()
        d = prov.to_dict()
        restored = OptimizationProvenance.from_dict(d)
        assert restored.organism == prov.organism
        assert restored.solver_backend == prov.solver_backend
        assert len(restored.decisions) == len(prov.decisions)
        assert restored.decisions[0].chosen_value == "ATG"
        assert restored.final_sequence == prov.final_sequence
        assert restored.solve_time_seconds == prov.solve_time_seconds
        assert restored.constraints_active == prov.constraints_active

    def test_from_dict_missing_keys_raises(self):
        with pytest.raises(ValueError, match="missing keys"):
            OptimizationProvenance.from_dict({"organism": "x"})

    def test_to_json_roundtrip(self):
        prov = self._make_provenance()
        json_str = prov.to_json()
        restored = OptimizationProvenance.from_json(json_str)
        assert restored.organism == prov.organism
        assert restored.solver_backend == prov.solver_backend

    def test_config_snapshot_preserved(self):
        prov = self._make_provenance()
        d = prov.to_dict()
        assert d["config_snapshot"]["gc_lo"] == 0.30
        assert d["config_snapshot"]["gc_hi"] == 0.70

    def test_constraints_active_list(self):
        prov = self._make_provenance()
        assert "GCInRange" in prov.constraints_active
        assert "NoCrypticSplice" in prov.constraints_active

    def test_repr(self):
        prov = self._make_provenance()
        r = repr(prov)
        assert "Homo_sapiens" in r
        assert "greedy" in r
        assert "1.230s" in r

    def test_decision_metadata_in_serialized_output(self):
        """Serialized provenance should include decision-level metadata."""
        prov = self._make_provenance()
        d = prov.to_dict()
        decision = d["decisions"][0]
        assert decision["decision_type"] == "codon_selected"
        assert decision["position"] == 0
        assert decision["chosen_value"] == "ATG"
        assert "rationale" in decision
        assert "constraint_context" in decision

    def test_empty_decisions(self):
        """OptimizationProvenance with no decisions is valid."""
        prov = OptimizationProvenance(
            input_protein="M",
            organism="E_coli",
            solver_backend="z3",
            config_snapshot={},
            decisions=[],
            final_sequence="ATG",
            solve_time_seconds=0.01,
            constraints_active=[],
        )
        assert len(prov.decisions) == 0
        d = prov.to_dict()
        assert d["decisions"] == []
