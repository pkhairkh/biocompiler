#!/usr/bin/env python3
"""
BioCompiler End-to-End Benchmark Suite
========================================
Comprehensive benchmark comparing BioCompiler vs DNAchisel across:
  1. DNA Optimization: CAI, GC content, restriction sites, speed
  2. Protein Folding/Structure: ESMFold pLDDT, solubility, stability
  3. All available gene sets × all supported organisms

This substitutes for wet-lab validation by using:
  - Publicly available gene sequences (UniProt/NCBI sourced)
  - Computational structure prediction (ESMFold heuristic)
  - CamSol solubility analysis
  - Statistical significance testing (paired t-test, Wilcoxon, Cohen's d)

Output:
  - CSV results file
  - JSON results file
  - Summary statistics
  - Per-organism breakdowns
"""

import csv
import json
import logging
import os
import sys
import time
from dataclasses import dataclass, field, asdict
from pathlib import Path
from typing import Any, Optional

# ── Ensure biocompiler is importable ──
REPO_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO_ROOT / "src"))

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(name)s: %(message)s",
)
logger = logging.getLogger("e2e_benchmark")

# ── Output directory ──
OUTPUT_DIR = Path("/home/z/my-project/download/e2e_benchmark")
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)


# ══════════════════════════════════════════════════════════════════════════════
# Data classes
# ══════════════════════════════════════════════════════════════════════════════

@dataclass
class DNAResult:
    """DNA optimization benchmark result for one gene × organism × tool."""
    gene_name: str
    organism: str
    target_organism: str
    protein_length: int
    tool: str  # "biocompiler" or "dnachisel"
    cai: float = 0.0
    gc_content: float = 0.0
    gc_std: float = 0.0
    restriction_site_count: int = 0
    cryptic_splice_sites: int = 0
    cpg_islands: int = 0
    mrna_stability: float = 0.0
    codon_pair_bias: float = 0.0
    runtime_s: float = 0.0
    success: bool = False
    error: str = ""
    gc_in_range: bool = False
    no_restriction_sites: bool = False
    cai_above_090: bool = False
    cai_above_095: bool = False


@dataclass
class FoldingResult:
    """Protein folding/structure benchmark result for one gene."""
    gene_name: str
    organism: str
    protein_length: int
    # ESMFold
    mean_plddt: float = 0.0
    plddt_method: str = ""
    plddt_confidence: str = ""
    # CamSol solubility
    solubility_score: float = 0.0
    solubility_class: str = ""
    intrinsic_score: float = 0.0
    aggregation_prone_regions: int = 0
    hydrophobicity_scale: str = ""
    # Stability
    mrna_stability_score: float = 0.0
    # Timing
    structure_time_s: float = 0.0
    solubility_time_s: float = 0.0
    success: bool = False
    error: str = ""


@dataclass
class BenchmarkRow:
    """Combined DNA + Folding result for CSV output."""
    gene_name: str
    source_organism: str
    target_organism: str
    protein_length: int
    # BioCompiler DNA
    bc_cai: float = 0.0
    bc_gc: float = 0.0
    bc_gc_std: float = 0.0
    bc_rs_count: int = 0
    bc_splice_sites: int = 0
    bc_cpg_islands: int = 0
    bc_mrna_stability: float = 0.0
    bc_runtime_s: float = 0.0
    bc_success: bool = False
    # DNAchisel DNA
    dc_cai: float = 0.0
    dc_gc: float = 0.0
    dc_gc_std: float = 0.0
    dc_rs_count: int = 0
    dc_splice_sites: int = 0
    dc_cpg_islands: int = 0
    dc_mrna_stability: float = 0.0
    dc_runtime_s: float = 0.0
    dc_success: bool = False
    # Folding
    mean_plddt: float = 0.0
    plddt_method: str = ""
    solubility_score: float = 0.0
    solubility_class: str = ""
    intrinsic_score: float = 0.0
    agg_regions: int = 0
    # Winner
    cai_winner: str = ""
    speed_winner: str = ""


# ══════════════════════════════════════════════════════════════════════════════
# Organism mapping
# ══════════════════════════════════════════════════════════════════════════════

# Map source organism names to supported target organisms for optimization
SOURCE_TO_TARGET: dict[str, list[str]] = {
    # Human → optimize for human, CHO, E. coli, yeast
    "Homo sapiens": ["Homo_sapiens", "CHO_K1", "Escherichia_coli", "Saccharomyces_cerevisiae"],
    # E. coli → optimize for E. coli, yeast
    "Escherichia coli": ["Escherichia_coli", "Saccharomyces_cerevisiae"],
    # Mouse → optimize for mouse, human, CHO
    "Mus musculus": ["Mus_musculus", "Homo_sapiens", "CHO_K1"],
    # Yeast → optimize for yeast, E. coli
    "Saccharomyces cerevisiae": ["Saccharomyces_cerevisiae", "Escherichia_coli"],
    # CHO → optimize for CHO, human
    "Cricetulus griseus": ["CHO_K1", "Homo_sapiens"],
    # Viral/pathogen genes → optimize for human, E. coli
    "Severe acute respiratory syndrome coronavirus 2": ["Homo_sapiens", "Escherichia_coli"],
    "Influenza A virus (A/Puerto Rico/8/1934 H1N1)": ["Homo_sapiens", "Escherichia_coli"],
    "Influenza A virus (A/Aichi/2/1968 H3N2)": ["Homo_sapiens", "Escherichia_coli"],
    "Human immunodeficiency virus type 1": ["Homo_sapiens", "Escherichia_coli"],
    "Human respiratory syncytial virus": ["Homo_sapiens", "Escherichia_coli"],
    "Rabies virus (strain PV)": ["Homo_sapiens", "Escherichia_coli"],
    "Zika virus": ["Homo_sapiens", "Escherichia_coli"],
    "Dengue virus type 2": ["Homo_sapiens", "Escherichia_coli"],
    "Mycobacterium tuberculosis": ["Escherichia_coli", "Homo_sapiens"],
    "Plasmodium falciparum": ["Homo_sapiens", "Escherichia_coli"],
    # Synthetic → all organisms
    "synthetic": ["Homo_sapiens", "Escherichia_coli", "Saccharomyces_cerevisiae", "CHO_K1", "Mus_musculus"],
}

# Normalise source organism name
_ORGANISM_ALIASES: dict[str, str] = {
    "homo sapiens": "Homo sapiens",
    "escherichia coli": "Escherichia coli",
    "mus musculus": "Mus musculus",
    "saccharomyces cerevisiae": "Saccharomyces cerevisiae",
    "cricetulus griseus": "Cricetulus griseus",
}

GC_RANGE = (0.30, 0.70)
ENZYME_PANEL = ["EcoRI", "BamHI", "HindIII", "XhoI", "NotI", "SalI", "PstI", "KpnI"]


def normalise_organism(name: str) -> str:
    """Normalise organism name with space format."""
    from biocompiler.organisms import SUPPORTED_ORGANISMS
    if name in SUPPORTED_ORGANISMS:
        return name.replace("_", " ")
    alias = _ORGANISM_ALIASES.get(name.lower())
    if alias:
        return alias
    return name.replace("_", " ")


def get_target_organisms(source_organism: str) -> list[str]:
    """Get list of target organisms for a given source organism."""
    # Try direct match
    targets = SOURCE_TO_TARGET.get(source_organism, [])
    if targets:
        return targets
    # Try normalised
    norm = normalise_organism(source_organism)
    targets = SOURCE_TO_TARGET.get(norm, [])
    if targets:
        return targets
    # Default: human + E. coli
    return ["Homo_sapiens", "Escherichia_coli"]


# ══════════════════════════════════════════════════════════════════════════════
# DNA optimization benchmark
# ══════════════════════════════════════════════════════════════════════════════

def run_biocompiler_dna(protein: str, organism: str, gene_name: str) -> DNAResult:
    """Run BioCompiler DNA optimization."""
    from biocompiler.optimization import optimize_sequence
    from biocompiler.benchmarking.metrics import (
        compute_cai_validated, compute_gc_distribution,
        count_restriction_sites, count_cryptic_splice_sites,
        count_cpg_islands, compute_mrna_stability_score,
    )

    result = DNAResult(
        gene_name=gene_name,
        organism="",
        target_organism=organism,
        protein_length=len(protein),
        tool="biocompiler",
    )

    t0 = time.perf_counter()
    try:
        opt = optimize_sequence(
            target_protein=protein,
            organism=organism,
            gc_lo=GC_RANGE[0],
            gc_hi=GC_RANGE[1],
        )
        seq = opt.sequence
        elapsed = time.perf_counter() - t0

        result.cai = compute_cai_validated(seq, organism)
        gc_profile = compute_gc_distribution(seq)
        result.gc_content = gc_profile.mean
        result.gc_std = gc_profile.std
        rs_dict = count_restriction_sites(seq, ENZYME_PANEL)
        result.restriction_site_count = sum(rs_dict.values())
        result.cryptic_splice_sites = count_cryptic_splice_sites(seq)
        result.cpg_islands = count_cpg_islands(seq)
        result.mrna_stability = compute_mrna_stability_score(seq)
        result.runtime_s = elapsed
        result.success = True
        result.gc_in_range = GC_RANGE[0] <= result.gc_content <= GC_RANGE[1]
        result.no_restriction_sites = result.restriction_site_count == 0
        result.cai_above_090 = result.cai >= 0.90
        result.cai_above_095 = result.cai >= 0.95

    except Exception as exc:
        result.runtime_s = time.perf_counter() - t0
        result.error = str(exc)[:200]
        logger.error("BioCompiler failed for %s/%s: %s", gene_name, organism, exc)

    return result


def run_dnachisel_dna(protein: str, organism: str, gene_name: str) -> DNAResult:
    """Run DNAchisel DNA optimization."""
    from biocompiler.benchmarking.dnachisel_adapter import DNAchiselAdapter, is_dnachisel_available
    from biocompiler.benchmarking.metrics import (
        compute_cai_validated, compute_gc_distribution,
        count_restriction_sites, count_cryptic_splice_sites,
        count_cpg_islands, compute_mrna_stability_score,
    )

    result = DNAResult(
        gene_name=gene_name,
        organism="",
        target_organism=organism,
        protein_length=len(protein),
        tool="dnachisel",
    )

    if not is_dnachisel_available():
        result.error = "DNAchisel not installed"
        return result

    t0 = time.perf_counter()
    try:
        adapter = DNAchiselAdapter()
        dc_result = adapter.optimize(
            protein=protein,
            organism=organism,
            constraints=[
                {"type": "gc_range", "gc_lo": GC_RANGE[0], "gc_hi": GC_RANGE[1]},
                {"type": "avoid_restriction", "enzymes": ENZYME_PANEL},
            ],
        )
        elapsed = time.perf_counter() - t0

        if dc_result.success:
            seq = dc_result.sequence
            result.cai = compute_cai_validated(seq, organism)
            gc_profile = compute_gc_distribution(seq)
            result.gc_content = gc_profile.mean
            result.gc_std = gc_profile.std
            rs_dict = count_restriction_sites(seq, ENZYME_PANEL)
            result.restriction_site_count = sum(rs_dict.values())
            result.cryptic_splice_sites = count_cryptic_splice_sites(seq)
            result.cpg_islands = count_cpg_islands(seq)
            result.mrna_stability = compute_mrna_stability_score(seq)
            result.runtime_s = elapsed
            result.success = True
            result.gc_in_range = GC_RANGE[0] <= result.gc_content <= GC_RANGE[1]
            result.no_restriction_sites = result.restriction_site_count == 0
            result.cai_above_090 = result.cai >= 0.90
            result.cai_above_095 = result.cai >= 0.95
        else:
            result.runtime_s = elapsed
            result.error = dc_result.error or "Unknown error"

    except Exception as exc:
        result.runtime_s = time.perf_counter() - t0
        result.error = str(exc)[:200]
        logger.error("DNAchisel failed for %s/%s: %s", gene_name, organism, exc)

    return result


# ══════════════════════════════════════════════════════════════════════════════
# Protein folding/structure benchmark
# ══════════════════════════════════════════════════════════════════════════════

def run_folding_analysis(protein: str, gene_name: str, organism: str = "Homo_sapiens") -> FoldingResult:
    """Run protein structure and solubility analysis."""
    from biocompiler.esmfold import predict_structure
    from biocompiler.camsol import compute_intrinsic_solubility

    result = FoldingResult(
        gene_name=gene_name,
        organism=organism,
        protein_length=len(protein),
    )

    # ESMFold structure prediction
    t0 = time.perf_counter()
    try:
        esm_result = predict_structure(protein, organism=organism)
        result.mean_plddt = esm_result.mean_plddt
        result.plddt_method = esm_result.method
        result.plddt_confidence = esm_result.confidence_level
        result.structure_time_s = time.perf_counter() - t0
    except Exception as exc:
        result.structure_time_s = time.perf_counter() - t0
        result.error = f"ESMFold: {str(exc)[:100]}"
        logger.warning("ESMFold failed for %s: %s", gene_name, exc)

    # CamSol solubility
    t1 = time.perf_counter()
    try:
        sol_result = compute_intrinsic_solubility(protein, organism=organism)
        result.solubility_score = sol_result.primary_score
        result.solubility_class = sol_result.classification
        result.intrinsic_score = sol_result.intrinsic_score
        result.aggregation_prone_regions = len(sol_result.aggregation_prone_regions)
        result.hydrophobicity_scale = sol_result.hydrophobicity_scale_used
        result.solubility_time_s = time.perf_counter() - t1
    except Exception as exc:
        result.solubility_time_s = time.perf_counter() - t1
        if result.error:
            result.error += f"; CamSol: {str(exc)[:80]}"
        else:
            result.error = f"CamSol: {str(exc)[:100]}"
        logger.warning("CamSol failed for %s: %s", gene_name, exc)

    # mRNA stability (from DNA-optimized result)
    try:
        from biocompiler.mrna_stability import score_mrna_stability
        result.mrna_stability_score = score_mrna_stability(protein, organism=organism)
    except Exception:
        pass

    result.success = True
    return result


# ══════════════════════════════════════════════════════════════════════════════
# Main benchmark runner
# ══════════════════════════════════════════════════════════════════════════════

def extract_protein(gene_data: dict) -> str:
    """Extract protein sequence from either legacy or new gene set format."""
    if "protein_sequence" in gene_data:
        return gene_data["protein_sequence"]
    if "protein" in gene_data:
        return gene_data["protein"]
    return ""


def extract_organism(gene_data: dict) -> str:
    """Extract organism from gene data."""
    return gene_data.get("organism", "synthetic")


def run_e2e_benchmark(
    gene_sets: dict[str, dict] | None = None,
    target_organisms_override: list[str] | None = None,
    skip_dnachisel: bool = False,
    skip_folding: bool = False,
    max_genes: int = 0,
) -> dict:
    """Run the complete end-to-end benchmark.

    Parameters
    ----------
    gene_sets : dict or None
        Gene sets to benchmark. If None, uses all available.
    target_organisms_override : list[str] or None
        If set, only test these target organisms.
    skip_dnachisel : bool
        Skip DNAchisel comparison.
    skip_folding : bool
        Skip folding/structure analysis.
    max_genes : int
        Max genes to process (0 = all).

    Returns
    -------
    dict with keys: "dna_results", "folding_results", "summary", "rows"
    """
    # ── Load gene sets ──
    if gene_sets is None:
        from biocompiler.benchmarking.gene_sets import get_all_gene_sets
        gene_sets = get_all_gene_sets()

    if max_genes > 0:
        items = list(gene_sets.items())[:max_genes]
        gene_sets = dict(items)

    total = len(gene_sets)
    logger.info("=" * 72)
    logger.info("BioCompiler E2E Benchmark")
    logger.info("=" * 72)
    logger.info("Total genes: %d", total)
    logger.info("DNAchisel: %s", "SKIP" if skip_dnachisel else "ENABLED")
    logger.info("Folding: %s", "SKIP" if skip_folding else "ENABLED")
    logger.info("Output: %s", OUTPUT_DIR)
    logger.info("=" * 72)

    dna_results: list[DNAResult] = []
    folding_results: list[FoldingResult] = []
    rows: list[BenchmarkRow] = []

    processed = 0
    for gene_name, gene_data in gene_sets.items():
        processed += 1
        protein = extract_protein(gene_data)
        source_org = extract_organism(gene_data)

        if not protein or len(protein) < 5:
            logger.warning("Skipping %s: empty/short protein", gene_name)
            continue

        # Determine target organisms
        if target_organisms_override:
            target_orgs = target_organisms_override
        else:
            target_orgs = get_target_organisms(source_org)

        # Skip viral genes for E. coli (membrane proteins, etc.)
        # that would fail with prokaryotic constraints

        for target_org in target_orgs:
            logger.info(
                "[%d/%d] %s (%d aa) → %s",
                processed, total, gene_name, len(protein), target_org,
            )

            # ── BioCompiler DNA ──
            bc_result = run_biocompiler_dna(protein, target_org, gene_name)
            dna_results.append(bc_result)

            # ── DNAchisel DNA ──
            dc_result: Optional[DNAResult] = None
            if not skip_dnachisel:
                dc_result = run_dnachisel_dna(protein, target_org, gene_name)
                dna_results.append(dc_result)

            # ── Folding (once per gene, not per target organism) ──
            fold_result: Optional[FoldingResult] = None
            if not skip_folding and target_org == target_orgs[0]:
                fold_result = run_folding_analysis(protein, gene_name, target_org)
                folding_results.append(fold_result)

            # ── Build row ──
            row = BenchmarkRow(
                gene_name=gene_name,
                source_organism=source_org,
                target_organism=target_org,
                protein_length=len(protein),
                bc_cai=bc_result.cai,
                bc_gc=bc_result.gc_content,
                bc_gc_std=bc_result.gc_std,
                bc_rs_count=bc_result.restriction_site_count,
                bc_splice_sites=bc_result.cryptic_splice_sites,
                bc_cpg_islands=bc_result.cpg_islands,
                bc_mrna_stability=bc_result.mrna_stability,
                bc_runtime_s=bc_result.runtime_s,
                bc_success=bc_result.success,
            )
            if dc_result:
                row.dc_cai = dc_result.cai
                row.dc_gc = dc_result.gc_content
                row.dc_gc_std = dc_result.gc_std
                row.dc_rs_count = dc_result.restriction_site_count
                row.dc_splice_sites = dc_result.cryptic_splice_sites
                row.dc_cpg_islands = dc_result.cpg_islands
                row.dc_mrna_stability = dc_result.mrna_stability
                row.dc_runtime_s = dc_result.runtime_s
                row.dc_success = dc_result.success
                # Winner determination
                if bc_result.success and dc_result.success:
                    row.cai_winner = "biocompiler" if bc_result.cai > dc_result.cai + 0.001 else (
                        "dnachisel" if dc_result.cai > bc_result.cai + 0.001 else "tie"
                    )
                    row.speed_winner = "biocompiler" if bc_result.runtime_s < dc_result.runtime_s * 0.95 else (
                        "dnachisel" if dc_result.runtime_s < bc_result.runtime_s * 0.95 else "tie"
                    )
            if fold_result:
                row.mean_plddt = fold_result.mean_plddt
                row.plddt_method = fold_result.plddt_method
                row.solubility_score = fold_result.solubility_score
                row.solubility_class = fold_result.solubility_class
                row.intrinsic_score = fold_result.intrinsic_score
                row.agg_regions = fold_result.aggregation_prone_regions

            rows.append(row)

    # ── Compute summary ──
    summary = compute_summary(dna_results, folding_results, rows)

    return {
        "dna_results": dna_results,
        "folding_results": folding_results,
        "rows": rows,
        "summary": summary,
    }


def compute_summary(
    dna_results: list[DNAResult],
    folding_results: list[FoldingResult],
    rows: list[BenchmarkRow],
) -> dict:
    """Compute aggregate statistics."""
    bc_dna = [r for r in dna_results if r.tool == "biocompiler" and r.success]
    dc_dna = [r for r in dna_results if r.tool == "dnachisel" and r.success]

    summary: dict[str, Any] = {}

    # DNA optimization counts
    summary["total_dna_tests"] = len(dna_results)
    summary["bc_success_count"] = len(bc_dna)
    summary["dc_success_count"] = len(dc_dna)
    summary["bc_success_rate"] = len(bc_dna) / max(1, len([r for r in dna_results if r.tool == "biocompiler"]))
    summary["dc_success_rate"] = len(dc_dna) / max(1, len([r for r in dna_results if r.tool == "dnachisel"]))

    # BioCompiler DNA aggregates
    if bc_dna:
        summary["bc_mean_cai"] = sum(r.cai for r in bc_dna) / len(bc_dna)
        summary["bc_median_cai"] = sorted(r.cai for r in bc_dna)[len(bc_dna) // 2]
        summary["bc_min_cai"] = min(r.cai for r in bc_dna)
        summary["bc_max_cai"] = max(r.cai for r in bc_dna)
        summary["bc_mean_gc"] = sum(r.gc_content for r in bc_dna) / len(bc_dna)
        summary["bc_mean_runtime"] = sum(r.runtime_s for r in bc_dna) / len(bc_dna)
        summary["bc_cai_above_090_rate"] = sum(1 for r in bc_dna if r.cai_above_090) / len(bc_dna)
        summary["bc_cai_above_095_rate"] = sum(1 for r in bc_dna if r.cai_above_095) / len(bc_dna)
        summary["bc_gc_in_range_rate"] = sum(1 for r in bc_dna if r.gc_in_range) / len(bc_dna)
        summary["bc_no_rs_rate"] = sum(1 for r in bc_dna if r.no_restriction_sites) / len(bc_dna)
    else:
        for k in ["bc_mean_cai", "bc_median_cai", "bc_min_cai", "bc_max_cai",
                   "bc_mean_gc", "bc_mean_runtime", "bc_cai_above_090_rate",
                   "bc_cai_above_095_rate", "bc_gc_in_range_rate", "bc_no_rs_rate"]:
            summary[k] = 0.0

    # DNAchisel DNA aggregates
    if dc_dna:
        summary["dc_mean_cai"] = sum(r.cai for r in dc_dna) / len(dc_dna)
        summary["dc_median_cai"] = sorted(r.cai for r in dc_dna)[len(dc_dna) // 2]
        summary["dc_min_cai"] = min(r.cai for r in dc_dna)
        summary["dc_max_cai"] = max(r.cai for r in dc_dna)
        summary["dc_mean_gc"] = sum(r.gc_content for r in dc_dna) / len(dc_dna)
        summary["dc_mean_runtime"] = sum(r.runtime_s for r in dc_dna) / len(dc_dna)
    else:
        for k in ["dc_mean_cai", "dc_median_cai", "dc_min_cai", "dc_max_cai",
                   "dc_mean_gc", "dc_mean_runtime"]:
            summary[k] = None

    # Head-to-head
    bc_wins_cai = sum(1 for r in rows if r.cai_winner == "biocompiler")
    dc_wins_cai = sum(1 for r in rows if r.cai_winner == "dnachisel")
    ties_cai = sum(1 for r in rows if r.cai_winner == "tie")
    summary["bc_cai_wins"] = bc_wins_cai
    summary["dc_cai_wins"] = dc_wins_cai
    summary["cai_ties"] = ties_cai

    bc_wins_speed = sum(1 for r in rows if r.speed_winner == "biocompiler")
    dc_wins_speed = sum(1 for r in rows if r.speed_winner == "dnachisel")
    ties_speed = sum(1 for r in rows if r.speed_winner == "tie")
    summary["bc_speed_wins"] = bc_wins_speed
    summary["dc_speed_wins"] = dc_wins_speed
    summary["speed_ties"] = ties_speed

    # Statistical significance (CAI)
    if bc_dna and dc_dna:
        try:
            from biocompiler.benchmarking.metrics import StatisticalComparison
            # Match by gene + organism
            bc_by_key = {}
            dc_by_key = {}
            for r in dna_results:
                key = (r.gene_name, r.target_organism)
                if r.tool == "biocompiler" and r.success:
                    bc_by_key[key] = r.cai
                elif r.tool == "dnachisel" and r.success:
                    dc_by_key[key] = r.cai
            matched_keys = set(bc_by_key.keys()) & set(dc_by_key.keys())
            if len(matched_keys) >= 2:
                bc_cais = [bc_by_key[k] for k in sorted(matched_keys)]
                dc_cais = [dc_by_key[k] for k in sorted(matched_keys)]
                comp = StatisticalComparison(bc_cais, dc_cais, "BioCompiler", "DNAchisel")
                stat = comp.summary()
                summary["cai_statistics"] = stat
                summary["cai_significant"] = stat["paired_t"]["significant"]
                summary["cai_cohens_d"] = stat["cohens_d"]
            else:
                summary["cai_statistics"] = None
                summary["cai_significant"] = False
                summary["cai_cohens_d"] = 0.0
        except Exception as exc:
            logger.warning("Statistical comparison failed: %s", exc)
            summary["cai_statistics"] = None
            summary["cai_significant"] = False
            summary["cai_cohens_d"] = 0.0

    # Per-organism breakdown
    org_breakdown: dict[str, dict] = {}
    for target_org in ["Homo_sapiens", "Escherichia_coli", "Saccharomyces_cerevisiae", "CHO_K1", "Mus_musculus"]:
        org_bc = [r for r in bc_dna if r.target_organism == target_org]
        org_dc = [r for r in dc_dna if r.target_organism == target_org]
        if org_bc:
            org_breakdown[target_org] = {
                "n_genes": len(org_bc),
                "bc_mean_cai": sum(r.cai for r in org_bc) / len(org_bc),
                "bc_mean_runtime": sum(r.runtime_s for r in org_bc) / len(org_bc),
                "dc_mean_cai": sum(r.cai for r in org_dc) / len(org_dc) if org_dc else None,
                "dc_mean_runtime": sum(r.runtime_s for r in org_dc) / len(org_dc) if org_dc else None,
            }
    summary["per_organism"] = org_breakdown

    # Folding aggregates
    if folding_results:
        success_folds = [f for f in folding_results if f.success]
        summary["total_folding_tests"] = len(folding_results)
        summary["folding_success_rate"] = len(success_folds) / len(folding_results)
        if success_folds:
            summary["mean_plddt"] = sum(f.mean_plddt for f in success_folds) / len(success_folds)
            summary["mean_solubility"] = sum(f.solubility_score for f in success_folds) / len(success_folds)
            summary["high_confidence_structure_rate"] = sum(
                1 for f in success_folds if f.plddt_confidence in ("high",)
            ) / len(success_folds)
            summary["soluble_rate"] = sum(
                1 for f in success_folds if f.solubility_class in ("highly_soluble", "soluble")
            ) / len(success_folds)
    else:
        summary["total_folding_tests"] = 0

    return summary


# ══════════════════════════════════════════════════════════════════════════════
# Report generation
# ══════════════════════════════════════════════════════════════════════════════

def save_csv(rows: list[BenchmarkRow], path: Path) -> None:
    """Save benchmark rows to CSV."""
    if not rows:
        return
    fieldnames = list(asdict(rows[0]).keys())
    with open(path, "w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        for row in rows:
            writer.writerow(asdict(row))
    logger.info("CSV results saved to %s", path)


def save_json(data: dict, path: Path) -> None:
    """Save benchmark data to JSON."""
    # Convert dataclasses to dicts
    def _serialize(obj):
        if hasattr(obj, '__dataclass_fields__'):
            return asdict(obj)
        if isinstance(obj, (list, tuple)):
            return [_serialize(x) for x in obj]
        if isinstance(obj, dict):
            return {k: _serialize(v) for k, v in obj.items()}
        if isinstance(obj, float):
            return round(obj, 6)
        return obj

    with open(path, "w", encoding="utf-8") as f:
        json.dump(_serialize(data), f, indent=2, default=str)
    logger.info("JSON results saved to %s", path)


def generate_report(summary: dict) -> str:
    """Generate a human-readable benchmark report."""
    lines = []
    lines.append("=" * 80)
    lines.append("  BioCompiler E2E Benchmark Report")
    lines.append("  DNA Optimization + Protein Folding/Structure Analysis")
    lines.append("=" * 80)
    lines.append("")

    # ── DNA Optimization ──
    lines.append("SECTION 1: DNA OPTIMIZATION (BioCompiler vs DNAchisel)")
    lines.append("-" * 80)
    lines.append(f"  Total DNA tests:        {summary.get('total_dna_tests', 0)}")
    lines.append(f"  BioCompiler successes:  {summary.get('bc_success_count', 0)} "
                 f"({summary.get('bc_success_rate', 0):.1%})")
    lines.append(f"  DNAchisel successes:    {summary.get('dc_success_count', 0) or 'N/A'} "
                 f"({summary.get('dc_success_rate', 0):.1%})" if summary.get('dc_success_rate') else
                 f"  DNAchisel successes:    N/A")
    lines.append("")

    # CAI comparison
    lines.append("  CAI (Codon Adaptation Index):")
    lines.append(f"    BioCompiler mean:     {summary.get('bc_mean_cai', 0):.4f}")
    lines.append(f"    BioCompiler median:   {summary.get('bc_median_cai', 0):.4f}")
    lines.append(f"    BioCompiler range:    [{summary.get('bc_min_cai', 0):.4f}, {summary.get('bc_max_cai', 0):.4f}]")
    if summary.get('dc_mean_cai') is not None:
        lines.append(f"    DNAchisel mean:       {summary['dc_mean_cai']:.4f}")
        lines.append(f"    DNAchisel median:     {summary.get('dc_median_cai', 0):.4f}")
        lines.append(f"    DNAchisel range:      [{summary.get('dc_min_cai', 0):.4f}, {summary.get('dc_max_cai', 0):.4f}]")
        cai_diff = summary['bc_mean_cai'] - summary['dc_mean_cai']
        lines.append(f"    Delta (BC - DC):      {cai_diff:+.4f}")
    lines.append("")

    # Quality rates
    lines.append("  Quality Rates (BioCompiler):")
    lines.append(f"    CAI >= 0.90:          {summary.get('bc_cai_above_090_rate', 0):.1%}")
    lines.append(f"    CAI >= 0.95:          {summary.get('bc_cai_above_095_rate', 0):.1%}")
    lines.append(f"    GC in range:          {summary.get('bc_gc_in_range_rate', 0):.1%}")
    lines.append(f"    No restriction sites: {summary.get('bc_no_rs_rate', 0):.1%}")
    lines.append("")

    # Speed
    lines.append("  Speed:")
    lines.append(f"    BioCompiler mean:     {summary.get('bc_mean_runtime', 0):.4f}s")
    if summary.get('dc_mean_runtime') is not None:
        lines.append(f"    DNAchisel mean:       {summary['dc_mean_runtime']:.4f}s")
        speed_ratio = summary['bc_mean_runtime'] / max(0.001, summary['dc_mean_runtime'])
        lines.append(f"    Ratio (BC/DC):        {speed_ratio:.2f}x")
    lines.append("")

    # Head-to-head
    lines.append("  Head-to-Head Results:")
    lines.append(f"    CAI:  BioCompiler wins {summary.get('bc_cai_wins', 0)}, "
                 f"DNAchisel wins {summary.get('dc_cai_wins', 0)}, "
                 f"ties {summary.get('cai_ties', 0)}")
    lines.append(f"    Speed: BioCompiler wins {summary.get('bc_speed_wins', 0)}, "
                 f"DNAchisel wins {summary.get('dc_speed_wins', 0)}, "
                 f"ties {summary.get('speed_ties', 0)}")
    lines.append("")

    # Statistical significance
    stats = summary.get("cai_statistics")
    if stats:
        lines.append("  Statistical Significance (CAI):")
        lines.append(f"    Paired t-test:  t={stats['paired_t']['statistic']:.3f}, "
                     f"p={stats['paired_t']['p_value']:.6f}, "
                     f"significant={'YES' if stats['paired_t']['significant'] else 'NO'}")
        lines.append(f"    Wilcoxon:       W={stats['wilcoxon']['statistic']:.3f}, "
                     f"p={stats['wilcoxon']['p_value']:.6f}, "
                     f"significant={'YES' if stats['wilcoxon']['significant'] else 'NO'}")
        lines.append(f"    Cohen's d:      {stats['cohens_d']:.3f}")
        lines.append(f"    Cliff's delta:  {stats['cliffs_delta']:.3f}")
        lines.append("")

    # Per-organism
    lines.append("  Per-Organism Breakdown:")
    for org, data in summary.get("per_organism", {}).items():
        lines.append(f"    {org}:")
        lines.append(f"      N genes: {data['n_genes']}")
        lines.append(f"      BC mean CAI:     {data['bc_mean_cai']:.4f}")
        lines.append(f"      BC mean runtime: {data['bc_mean_runtime']:.4f}s")
        if data.get('dc_mean_cai') is not None:
            lines.append(f"      DC mean CAI:     {data['dc_mean_cai']:.4f}")
            lines.append(f"      DC mean runtime: {data['dc_mean_runtime']:.4f}s")
        lines.append("")

    # ── Protein Folding ──
    lines.append("SECTION 2: PROTEIN FOLDING & STRUCTURE ANALYSIS")
    lines.append("-" * 80)
    total_fold = summary.get("total_folding_tests", 0)
    lines.append(f"  Total folding tests:    {total_fold}")
    if total_fold > 0:
        lines.append(f"  Success rate:           {summary.get('folding_success_rate', 0):.1%}")
        lines.append(f"  Mean pLDDT:             {summary.get('mean_plddt', 0):.1f}")
        lines.append(f"  Mean solubility:        {summary.get('mean_solubility', 0):.4f}")
        lines.append(f"  High-confidence struct: {summary.get('high_confidence_structure_rate', 0):.1%}")
        lines.append(f"  Soluble classification: {summary.get('soluble_rate', 0):.1%}")
    lines.append("")

    lines.append("=" * 80)
    return "\n".join(lines)


# ══════════════════════════════════════════════════════════════════════════════
# Entry point
# ══════════════════════════════════════════════════════════════════════════════

def main():
    """Run the E2E benchmark and save results."""
    logger.info("Starting BioCompiler E2E Benchmark...")

    results = run_e2e_benchmark(
        skip_dnachisel=False,
        skip_folding=False,
        max_genes=0,  # 0 = all genes
    )

    # Save CSV
    csv_path = OUTPUT_DIR / "e2e_benchmark_results.csv"
    save_csv(results["rows"], csv_path)

    # Save JSON
    json_path = OUTPUT_DIR / "e2e_benchmark_results.json"
    save_json(results, json_path)

    # Generate and save report
    report = generate_report(results["summary"])
    report_path = OUTPUT_DIR / "e2e_benchmark_report.txt"
    with open(report_path, "w") as f:
        f.write(report)

    # Print report
    print(report)

    # Also save DNA-only and folding-only JSONs
    save_json(
        {"results": results["dna_results"]},
        OUTPUT_DIR / "dna_optimization_results.json",
    )
    save_json(
        {"results": results["folding_results"]},
        OUTPUT_DIR / "folding_results.json",
    )

    logger.info("Benchmark complete. Results in %s", OUTPUT_DIR)
    return results


if __name__ == "__main__":
    main()
