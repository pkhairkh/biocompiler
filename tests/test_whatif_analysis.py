"""Tests for WhatIfScenario, WhatIfAnalyzer, and WhatIfReport.

Covers:
- WhatIfScenario creation and validation
- WhatIfScenario serialization/deserialization
- WhatIfAnalyzer initialization and validation
- GC relaxation/tightening analysis
- Constraint removal analysis
- Organism switch analysis
- WhatIfReport generation
- JSON export/import
"""

import json
import os
import tempfile
import pytest
from biocompiler.whatif_analysis import WhatIfScenario, WhatIfAnalyzer, WhatIfReport


class TestWhatIfScenario:
    def test_basic_creation(self):
        scenario = WhatIfScenario(
            description="Test scenario",
            parameter_changed="gc_hi",
            original_value=0.70,
            alternative_value=0.75,
            predicted_cai=0.85,
            predicted_gc=0.65,
            constraint_satisfaction={"gc_in_range": True},
            feasibility="feasible",
            baseline_cai=0.80,
            baseline_gc=0.60,
        )
        assert scenario.description == "Test scenario"
        assert scenario.parameter_changed == "gc_hi"
        assert scenario.feasibility == "feasible"
        assert scenario.predicted_cai == pytest.approx(0.85, rel=1e-6)

    def test_cai_delta_computed(self):
        scenario = WhatIfScenario(
            description="Test",
            parameter_changed="gc_hi",
            original_value=0.70,
            alternative_value=0.75,
            predicted_cai=0.85,
            predicted_gc=0.65,
            constraint_satisfaction=None,
            feasibility="feasible",
            baseline_cai=0.80,
        )
        assert scenario.cai_delta == pytest.approx(0.05, abs=1e-6)

    def test_timestamp_auto_set(self):
        scenario = WhatIfScenario(
            description="Test",
            parameter_changed="gc_hi",
            original_value=0.70,
            alternative_value=0.75,
            predicted_cai=None,
            predicted_gc=None,
            constraint_satisfaction=None,
            feasibility="unknown",
        )
        assert scenario.timestamp != ""

    def test_invalid_feasibility_raises(self):
        with pytest.raises(ValueError, match="feasibility"):
            WhatIfScenario(
                description="Test",
                parameter_changed="gc_hi",
                original_value=0.70,
                alternative_value=0.75,
                predicted_cai=None,
                predicted_gc=None,
                constraint_satisfaction=None,
                feasibility="invalid_value",
            )

    def test_to_dict(self):
        scenario = WhatIfScenario(
            description="Test",
            parameter_changed="gc_hi",
            original_value=0.70,
            alternative_value=0.75,
            predicted_cai=0.85,
            predicted_gc=0.65,
            constraint_satisfaction={"gc_in_range": True},
            feasibility="feasible",
        )
        d = scenario.to_dict()
        assert isinstance(d, dict)
        assert d["description"] == "Test"
        assert d["predicted_cai"] == pytest.approx(0.85, rel=1e-6)
        assert d["feasibility"] == "feasible"

    def test_from_dict(self):
        d = {
            "description": "Test",
            "parameter_changed": "gc_hi",
            "original_value": 0.70,
            "alternative_value": 0.75,
            "predicted_cai": 0.85,
            "predicted_gc": 0.65,
            "constraint_satisfaction": {"gc_in_range": True},
            "feasibility": "feasible",
            "baseline_cai": 0.80,
        }
        scenario = WhatIfScenario.from_dict(d)
        assert scenario.description == "Test"
        assert scenario.predicted_cai == pytest.approx(0.85, rel=1e-6)

    def test_from_dict_missing_keys_raises(self):
        d = {"description": "Test"}  # Missing required keys
        with pytest.raises(ValueError, match="missing keys"):
            WhatIfScenario.from_dict(d)

    def test_repr(self):
        scenario = WhatIfScenario(
            description="Test",
            parameter_changed="gc_hi",
            original_value=0.70,
            alternative_value=0.75,
            predicted_cai=0.85,
            predicted_gc=0.65,
            constraint_satisfaction=None,
            feasibility="feasible",
        )
        r = repr(scenario)
        assert "WhatIfScenario" in r
        assert "feasible" in r

    def test_none_predicted_cai(self):
        scenario = WhatIfScenario(
            description="Test",
            parameter_changed="gc_hi",
            original_value=0.70,
            alternative_value=0.75,
            predicted_cai=None,
            predicted_gc=None,
            constraint_satisfaction=None,
            feasibility="unknown",
        )
        assert scenario.predicted_cai is None
        assert scenario.cai_delta is None


class TestWhatIfAnalyzer:
    def test_init_human(self):
        analyzer = WhatIfAnalyzer(protein="MVSKGE", organism="Homo_sapiens")
        assert analyzer.organism == "Homo_sapiens"
        assert analyzer.protein == "MVSKGE"

    def test_init_ecoli(self):
        analyzer = WhatIfAnalyzer(protein="MVSKGE", organism="Escherichia_coli")
        assert analyzer.organism == "Escherichia_coli"

    def test_unsupported_organism_raises(self):
        with pytest.raises(ValueError, match="Unsupported organism"):
            WhatIfAnalyzer(protein="MVSKGE", organism="Alien_martian")

    def test_gc_lo_hi_stored(self):
        analyzer = WhatIfAnalyzer(protein="MVSKGE", organism="Homo_sapiens", gc_lo=0.35, gc_hi=0.65)
        assert analyzer.gc_lo == pytest.approx(0.35, rel=1e-6)
        assert analyzer.gc_hi == pytest.approx(0.65, rel=1e-6)

    def test_analyze_gc_relaxation(self):
        analyzer = WhatIfAnalyzer(protein="MVSKGE", organism="Escherichia_coli")
        dna = "ATGGTTTCTAAAGGTGAA"
        scenario = analyzer.analyze_gc_relaxation(dna, current_gc_hi=0.70, alternative_gc_hi=0.80)
        assert isinstance(scenario, WhatIfScenario)
        assert scenario.parameter_changed == "gc_hi"
        assert scenario.feasibility in ("feasible", "infeasible", "unknown")

    def test_analyze_gc_tightening(self):
        analyzer = WhatIfAnalyzer(protein="MVSKGE", organism="Escherichia_coli")
        dna = "ATGGTTTCTAAAGGTGAA"
        scenario = analyzer.analyze_gc_tightening(dna, 0.30, 0.70, 0.35, 0.65)
        assert isinstance(scenario, WhatIfScenario)
        assert scenario.parameter_changed == "gc_range"

    def test_analyze_constraint_removal_restriction_sites(self):
        analyzer = WhatIfAnalyzer(protein="MVSKGE", organism="Escherichia_coli")
        dna = "ATGGTTTCTAAAGGTGAA"
        scenario = analyzer.analyze_constraint_removal(dna, "restriction_sites")
        assert isinstance(scenario, WhatIfScenario)
        assert scenario.parameter_changed == "restriction_sites"

    def test_analyze_organism_switch(self):
        analyzer = WhatIfAnalyzer(protein="MVSKGE", organism="Escherichia_coli")
        dna = "ATGGTTTCTAAAGGTGAA"
        scenario = analyzer.analyze_organism_switch(dna, "Homo_sapiens")
        assert isinstance(scenario, WhatIfScenario)
        assert scenario.parameter_changed == "organism"

    def test_analyze_organism_switch_unsupported_raises(self):
        analyzer = WhatIfAnalyzer(protein="MVSKGE", organism="Escherichia_coli")
        dna = "ATGGTTTCTAAAGGTGAA"
        with pytest.raises(ValueError, match="Unsupported organism"):
            analyzer.analyze_organism_switch(dna, "Alien_martian")


class TestWhatIfReport:
    def test_generate_empty_scenarios(self):
        report = WhatIfReport.generate([])
        assert "No scenarios" in report

    def test_generate_with_scenarios(self):
        scenarios = [
            WhatIfScenario(
                description="Test scenario 1",
                parameter_changed="gc_hi",
                original_value=0.70,
                alternative_value=0.75,
                predicted_cai=0.85,
                predicted_gc=0.65,
                constraint_satisfaction={"gc_in_range": True},
                feasibility="feasible",
                baseline_cai=0.80,
            ),
        ]
        report = WhatIfReport.generate(scenarios)
        assert "What-If Analysis Report" in report
        assert "gc_hi" in report

    def test_to_json(self):
        scenarios = [
            WhatIfScenario(
                description="Test",
                parameter_changed="gc_hi",
                original_value=0.70,
                alternative_value=0.75,
                predicted_cai=0.85,
                predicted_gc=0.65,
                constraint_satisfaction=None,
                feasibility="feasible",
            ),
        ]
        with tempfile.NamedTemporaryFile(suffix=".json", delete=False) as f:
            filepath = f.name
        try:
            WhatIfReport.to_json(scenarios, filepath)
            assert os.path.exists(filepath)
            with open(filepath) as f:
                data = json.load(f)
            assert data["report_type"] == "whatif_analysis"
            assert len(data["scenarios"]) == 1
        finally:
            os.unlink(filepath)

    def test_from_json(self):
        scenarios = [
            WhatIfScenario(
                description="Test",
                parameter_changed="gc_hi",
                original_value=0.70,
                alternative_value=0.75,
                predicted_cai=0.85,
                predicted_gc=0.65,
                constraint_satisfaction=None,
                feasibility="feasible",
            ),
        ]
        with tempfile.NamedTemporaryFile(suffix=".json", delete=False) as f:
            filepath = f.name
        try:
            WhatIfReport.to_json(scenarios, filepath)
            loaded = WhatIfReport.from_json(filepath)
            assert len(loaded) == 1
            assert loaded[0].description == "Test"
        finally:
            os.unlink(filepath)

    def test_from_json_invalid_file_raises(self):
        with tempfile.NamedTemporaryFile(suffix=".json", delete=False, mode='w') as f:
            f.write('{"not_scenarios": []}')
            filepath = f.name
        try:
            with pytest.raises(ValueError, match="missing 'scenarios'"):
                WhatIfReport.from_json(filepath)
        finally:
            os.unlink(filepath)
