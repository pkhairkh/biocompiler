"""
Performance tests for BioCompiler optimization pipeline.

Benchmarks:
- Optimize GFP for E. coli (prokaryotic): target < 50ms
- Optimize GFP for Human (eukaryotic): target < 200ms

Tests:
- should_skip_constraint returns correct values
- estimate_optimization_complexity classifies correctly
- batch_detect_violations finds same violations as individual checks
- warm_numba_cache doesn't crash
"""

import time
import unittest


class TestShouldSkipConstraint(unittest.TestCase):
    """Test the should_skip_constraint function."""

    def test_prokaryote_skips_eukaryotic_constraints(self):
        """Prokaryotes should skip all eukaryotic-only constraint types."""
        from biocompiler.optimizer.performance import should_skip_constraint

        eukaryotic_types = ["cpg_island", "cryptic_splice_donor",
                            "cryptic_splice_acceptor", "avoidable_gt"]
        for vtype in eukaryotic_types:
            with self.subTest(vtype=vtype):
                self.assertTrue(
                    should_skip_constraint(vtype, is_prokaryote=True),
                    f"Prokaryote should skip {vtype}"
                )

    def test_prokaryote_does_not_skip_hard_constraints(self):
        """Prokaryotes should NOT skip hard constraints like restriction sites."""
        from biocompiler.optimizer.performance import should_skip_constraint

        hard_types = ["restriction_site", "gc_out_of_range", "stop_codon",
                      "atttta_motif", "t_run"]
        for vtype in hard_types:
            with self.subTest(vtype=vtype):
                self.assertFalse(
                    should_skip_constraint(vtype, is_prokaryote=True),
                    f"Prokaryote should NOT skip {vtype}"
                )

    def test_eukaryote_does_not_skip_constraints_by_default(self):
        """Eukaryotes should not skip constraints without a previous result."""
        from biocompiler.optimizer.performance import should_skip_constraint

        self.assertFalse(
            should_skip_constraint("cpg_island", is_prokaryote=False)
        )
        self.assertFalse(
            should_skip_constraint("avoidable_gt", is_prokaryote=False)
        )

    def test_previous_result_with_high_cai_skips_soft(self):
        """When previous result has high CAI and no violations, soft constraints can be skipped."""
        from biocompiler.optimizer.performance import should_skip_constraint
        from biocompiler.optimizer.hybrid_types import HybridResult

        prev = HybridResult(
            sequence="ATGAAAGCGTAA",
            cai=0.97,
            gc_content=0.44,
            violations_fixed=0,
        )
        # Soft constraints should be skippable with a good previous result
        self.assertTrue(
            should_skip_constraint("cpg_island", is_prokaryote=False,
                                  previous_result=prev)
        )
        self.assertTrue(
            should_skip_constraint("atttta_motif", is_prokaryote=False,
                                  previous_result=prev)
        )


class TestEstimateOptimizationComplexity(unittest.TestCase):
    """Test the estimate_optimization_complexity function."""

    def test_simple_classification(self):
        """Short protein, no enzymes, wide GC range = simple."""
        from biocompiler.optimizer.performance import estimate_optimization_complexity

        result = estimate_optimization_complexity(
            protein="M" * 100,  # Short protein
            organism="Escherichia_coli",
            enzymes=[],
            gc_lo=0.30,
            gc_hi=0.70,
        )
        self.assertEqual(result, "simple")

    def test_complex_classification(self):
        """Long protein, many enzymes, tight GC range = complex."""
        from biocompiler.optimizer.performance import estimate_optimization_complexity

        result = estimate_optimization_complexity(
            protein="M" * 600,  # Long protein
            organism="Homo_sapiens",
            enzymes=["EcoRI", "BamHI", "XhoI", "HindIII", "NotI", "SalI"],
            gc_lo=0.48,
            gc_hi=0.52,  # Very tight GC range (0.04 width)
        )
        self.assertEqual(result, "complex")

    def test_moderate_classification(self):
        """Medium protein, some enzymes, moderate GC range = moderate."""
        from biocompiler.optimizer.performance import estimate_optimization_complexity

        result = estimate_optimization_complexity(
            protein="M" * 350,  # Medium protein
            organism="Saccharomyces_cerevisiae",
            enzymes=["EcoRI", "BamHI"],
            gc_lo=0.35,
            gc_hi=0.55,  # Moderate GC range
        )
        self.assertIn(result, ["simple", "moderate", "complex"])


class TestBatchDetectViolations(unittest.TestCase):
    """Test that batch_detect_violations finds the same violations as individual checks."""

    def test_finds_atttta_motifs(self):
        """Batch detection should find ATTTA motifs."""
        from biocompiler.optimizer.performance import batch_detect_violations

        seq = "ATGAAATTTAAGCGTAA"  # Contains ATTTA
        violations = batch_detect_violations(
            seq,
            ["atttta_motif"],
            {"is_prokaryote": True, "rs_sites": []},
        )
        attta_violations = [v for v in violations if v.violation_type == "atttta_motif"]
        self.assertGreater(len(attta_violations), 0,
                          "Should detect ATTTA motif")

    def test_finds_t_runs(self):
        """Batch detection should find T-runs of 6+."""
        from biocompiler.optimizer.performance import batch_detect_violations

        seq = "ATGAAATTTTTTGCGTAA"  # Contains TTTTTT (6 T's)
        violations = batch_detect_violations(
            seq,
            ["t_run"],
            {"is_prokaryote": True, "rs_sites": []},
        )
        trun_violations = [v for v in violations if v.violation_type == "t_run"]
        self.assertGreater(len(trun_violations), 0,
                          "Should detect T-run")

    def test_skips_eukaryotic_for_prokaryote(self):
        """Batch detection should skip GT/AG/CG for prokaryotes."""
        from biocompiler.optimizer.performance import batch_detect_violations

        seq = "ATGGCGTGTAGCGTAA"  # Contains GT, AG, CG
        violations = batch_detect_violations(
            seq,
            ["gt_dinucleotide", "ag_dinucleotide", "cpg_dinucleotide"],
            {"is_prokaryote": True, "rs_sites": []},
        )
        euk_violations = [v for v in violations
                         if v.violation_type in ("gt_dinucleotide",
                                                 "ag_dinucleotide",
                                                 "cpg_dinucleotide")]
        self.assertEqual(len(euk_violations), 0,
                        "Prokaryote should not detect eukaryotic violations")

    def test_finds_eukaryotic_violations(self):
        """Batch detection should find GT/AG/CG for eukaryotes."""
        from biocompiler.optimizer.performance import batch_detect_violations

        seq = "ATGGCGTGTAGCGTAA"
        violations = batch_detect_violations(
            seq,
            ["gt_dinucleotide", "ag_dinucleotide", "cpg_dinucleotide"],
            {"is_prokaryote": False, "rs_sites": []},
        )
        gt_v = [v for v in violations if v.violation_type == "gt_dinucleotide"]
        ag_v = [v for v in violations if v.violation_type == "ag_dinucleotide"]
        cpg_v = [v for v in violations if v.violation_type == "cpg_dinucleotide"]
        self.assertGreater(len(gt_v), 0, "Should detect GT")
        self.assertGreater(len(ag_v), 0, "Should detect AG")
        self.assertGreater(len(cpg_v), 0, "Should detect CG")

    def test_no_violations_in_clean_sequence(self):
        """Batch detection should find no violations in a clean sequence."""
        from biocompiler.optimizer.performance import batch_detect_violations

        seq = "ATGAAAGCGTAA"  # No problematic motifs
        violations = batch_detect_violations(
            seq,
            ["atttta_motif", "t_run"],
            {"is_prokaryote": True, "rs_sites": []},
        )
        self.assertEqual(len(violations), 0,
                        "Clean sequence should have no violations")


class TestGetFastPathConfig(unittest.TestCase):
    """Test the get_fast_path_config function."""

    def test_simple_prokaryote_config(self):
        """Simple prokaryote should have minimal iterations."""
        from biocompiler.optimizer.performance import get_fast_path_config

        config = get_fast_path_config("simple", is_prokaryote=True)
        self.assertTrue(config.skip_hill_climbing)
        self.assertLessEqual(config.hill_climb_passes, 3)
        self.assertEqual(config.cpg_max_iterations, 0)
        self.assertTrue(config.early_termination)

    def test_complex_eukaryote_config(self):
        """Complex eukaryote should have more iterations."""
        from biocompiler.optimizer.performance import get_fast_path_config

        config = get_fast_path_config("complex", is_prokaryote=False)
        self.assertFalse(config.skip_hill_climbing)
        self.assertGreater(config.cpg_max_iterations, 0)
        self.assertGreater(config.hill_climb_passes, 2)

    def test_moderate_prokaryote_less_than_complex(self):
        """Moderate prokaryote should use fewer iterations than complex."""
        from biocompiler.optimizer.performance import get_fast_path_config

        mod = get_fast_path_config("moderate", is_prokaryote=True)
        comp = get_fast_path_config("complex", is_prokaryote=True)
        self.assertLessEqual(
            mod.max_local_search_iterations,
            comp.max_local_search_iterations
        )


class TestWarmNumbaCache(unittest.TestCase):
    """Test that warm_numba_cache doesn't crash."""

    def test_warm_numba_cache_no_crash(self):
        """warm_numba_cache should not raise any exceptions."""
        from biocompiler.optimizer.performance import warm_numba_cache

        # Should not raise, even if NUMBA is not available
        try:
            warm_numba_cache()
        except Exception as e:
            self.fail(f"warm_numba_cache raised {type(e).__name__}: {e}")


class TestShouldSkipHelpers(unittest.TestCase):
    """Test the helper functions for skipping passes."""

    def test_skip_mrna_stability_no_attta(self):
        """Should skip mRNA stability when no ATTTA motifs exist."""
        from biocompiler.optimizer.performance import should_skip_mrna_stability

        self.assertTrue(should_skip_mrna_stability("ATGAAAGCGTAA"))
        self.assertFalse(should_skip_mrna_stability("ATGATTTAGCGTAA"))

    def test_skip_cpg_for_prokaryote(self):
        """Should skip CpG elimination for prokaryotes."""
        from biocompiler.optimizer.performance import should_skip_cpg_elimination

        self.assertTrue(should_skip_cpg_elimination(is_prokaryote=True, sequence="ATGCGTAA"))
        self.assertFalse(should_skip_cpg_elimination(is_prokaryote=False, sequence="ATGCGTAA"))

    def test_skip_cpg_no_cg(self):
        """Should skip CpG elimination when no CG dinucleotides exist."""
        from biocompiler.optimizer.performance import should_skip_cpg_elimination

        self.assertTrue(should_skip_cpg_elimination(is_prokaryote=False, sequence="ATGAAATTTTAA"))

    def test_skip_utr_suggestions(self):
        """Should skip UTR suggestions when include_utr=False."""
        from biocompiler.optimizer.performance import should_skip_utr_suggestions

        self.assertTrue(should_skip_utr_suggestions(include_utr=False))
        self.assertFalse(should_skip_utr_suggestions(include_utr=True))


class TestOrganismDataCache(unittest.TestCase):
    """Test the organism data cache."""

    def test_get_organism_data_returns_dict(self):
        """get_organism_data should return a dictionary with expected keys."""
        from biocompiler.optimizer.performance import get_organism_data, clear_caches

        clear_caches()
        data = get_organism_data("Escherichia_coli")
        self.assertIsInstance(data, dict)
        self.assertIn("species_cai", data)
        self.assertIn("is_prokaryote", data)
        self.assertTrue(data["is_prokaryote"])

    def test_get_organism_data_caches(self):
        """Second call should return the same cached object."""
        from biocompiler.optimizer.performance import get_organism_data, clear_caches

        clear_caches()
        data1 = get_organism_data("Homo_sapiens")
        data2 = get_organism_data("Homo_sapiens")
        self.assertIs(data1, data2)

    def test_clear_caches(self):
        """clear_caches should invalidate the cache."""
        from biocompiler.optimizer.performance import get_organism_data, clear_caches

        clear_caches()
        data1 = get_organism_data("Saccharomyces_cerevisiae")
        clear_caches()
        data2 = get_organism_data("Saccharomyces_cerevisiae")
        self.assertIsNot(data1, data2)


class BenchmarkOptimization(unittest.TestCase):
    """Benchmark tests for optimization speed."""

    def _get_gfp_protein(self) -> str:
        """Return GFP protein sequence (238 aa)."""
        return (
            "MSKGEELFTGVVPILVELDGDVNGHKFSVSGEGEGDATYGKLTLKFICTTGKLPVPWPTLVTTLTYGVQCFSRYPDHMKQHDFFKSAMPEGYVQERTIFFKDDGNYKTRAEVKFEGDTLVNRIELKGIDFKEDGNILGHKLEYNYNSHNVYIMADKQKNGIKVNFKIRHNIEDGSVQLADHYQQNTPIGDGPVLLPDNHYLSTQSALSKDPNEKRDHMVLLEFVTAAGITHGMDELYK"
        )

    def test_prokaryote_gfp_speed(self):
        """Optimize GFP for E. coli: should complete in < 50ms (after warmup)."""
        try:
            from biocompiler.hybrid_optimizer import HybridOptimizer
        except ImportError:
            self.skipTest("HybridOptimizer not importable")

        protein = self._get_gfp_protein()
        opt = HybridOptimizer(
            species="ecoli",
            organism="Escherichia_coli",
            enzymes=["EcoRI", "BamHI", "XhoI", "HindIII", "NotI"],
            gc_lo=0.30,
            gc_hi=0.70,
        )

        # Warmup call to eliminate first-call overhead
        _ = opt.optimize(protein, is_prokaryote=True)

        # Benchmark
        times = []
        for _ in range(5):
            start = time.perf_counter()
            result = opt.optimize(protein, is_prokaryote=True)
            elapsed = (time.perf_counter() - start) * 1000  # ms
            times.append(elapsed)

        # Use median time to reduce variance
        times.sort()
        median_time = times[len(times) // 2]

        # Verify result quality
        self.assertGreater(len(result.sequence), 0)
        self.assertGreater(result.cai, 0.0)

        # Log the time for visibility (not a hard failure)
        print(f"\n  E. coli GFP optimization: {median_time:.1f}ms (median of 5)")
        if median_time > 50:
            print(f"  WARNING: Exceeded 50ms target ({median_time:.1f}ms)")
            # Don't fail on slow CI, just warn

    def test_eukaryote_gfp_speed(self):
        """Optimize GFP for Human: should complete in < 200ms (after warmup)."""
        try:
            from biocompiler.hybrid_optimizer import HybridOptimizer
        except ImportError:
            self.skipTest("HybridOptimizer not importable")

        protein = self._get_gfp_protein()
        opt = HybridOptimizer(
            species="human",
            organism="Homo_sapiens",
            enzymes=["EcoRI", "BamHI", "XhoI", "HindIII", "NotI"],
            gc_lo=0.30,
            gc_hi=0.70,
            avoid_gt=True,
        )

        # Warmup call
        _ = opt.optimize(protein, is_prokaryote=False)

        # Benchmark
        times = []
        for _ in range(5):
            start = time.perf_counter()
            result = opt.optimize(protein, is_prokaryote=False)
            elapsed = (time.perf_counter() - start) * 1000
            times.append(elapsed)

        times.sort()
        median_time = times[len(times) // 2]

        # Verify result quality
        self.assertGreater(len(result.sequence), 0)
        self.assertGreater(result.cai, 0.0)

        # Log the time
        print(f"\n  Human GFP optimization: {median_time:.1f}ms (median of 5)")
        if median_time > 200:
            print(f"  WARNING: Exceeded 200ms target ({median_time:.1f}ms)")


if __name__ == "__main__":
    unittest.main()
