#!/usr/bin/env python3
"""
Mega Benchmark: BioCompiler vs DNAchisel — Full Gene Panel Validation
======================================================================

Comprehensive benchmark across ALL available gene sets:
  - E. coli Extended (20 genes)
  - Human Therapeutic (10 genes)
  - Human Signaling (4 genes)
  - Yeast Industrial (4 genes)
  - Vaccine Antigens (12 genes)
  - Stress Test (5 genes)
  - Published Validation Sequences (with known CAI values)

Total: 55+ genes across 3+ organisms

Validation against public data:
  - Sharp & Li (1987) published CAI values
  - Puigbo et al. (2008) CAIcal server values
  - NCBI RefSeq coding sequences as ground truth
  - Expression-level rank correlation (high/medium/low CAI ordering)

Output: JSON + CSV + PNG plots + PDF-ready summary
"""

from __future__ import annotations

import json
import math
import os
import sys
import time
from dataclasses import dataclass, field, asdict
from pathlib import Path
from typing import Any

# ── Project imports ──
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "src"))

from biocompiler.optimizer import optimize_sequence
from biocompiler.benchmarking.metrics import (
    compute_cai_validated,
    compute_gc_distribution,
    count_restriction_sites,
    count_cryptic_splice_sites,
    compute_mrna_stability_score,
    count_cpg_islands,
    STANDARD_ENZYME_PANEL,
)
from biocompiler.benchmarking.gene_sets import (
    E_COLI_EXTENDED,
    HUMAN_THERAPEUTIC_GENES,
    HUMAN_THERAPEUTIC,
    HUMAN_SIGNALING,
    YEAST_INDUSTRIAL,
    VACCINE_ANTIGEN_GENES,
    STRESS_TEST_GENES,
)
from biocompiler.benchmarking.cai_published_values import (
    PUBLISHED_CAI_VALUES,
    VALIDATION_SEQUENCES,
)
from biocompiler.sequence.scanner import gc_content

# ── DNAchisel imports ──
_DNACHISEL_AVAILABLE = False
try:
    from dnachisel import (
        DnaOptimizationProblem,
        AvoidPattern,
        CodonOptimize,
        EnforceGCContent,
        EnforceTranslation,
    )
    _DNACHISEL_AVAILABLE = True
except ImportError:
    pass

from biocompiler.organisms import PREFERRED_CODON_TABLES, CODON_ADAPTIVENESS_TABLES
from biocompiler.shared.constants import AA_TO_CODONS

_SPECIES_MAP = {
    "Escherichia_coli": "e_coli",
    "Homo_sapiens": "h_sapiens",
    "Saccharomyces_cerevisiae": "s_cerevisiae",
    "Mus_musculus": "m_musculus",
}

GC_LO = 0.30
GC_HI = 0.70


# ═══════════════════════════════════════════════════════════════════════
# Data Classes
# ═══════════════════════════════════════════════════════════════════════

@dataclass
class GeneResult:
    gene_name: str
    organism: str
    tool: str
    protein_length: int
    cai: float
    gc_content: float
    gc_std: float
    runtime_ms: float
    success: bool
    error: str = ""
    restriction_site_total: int = 0
    cryptic_splice_sites: int = 0
    mrna_stability: float = 0.0
    cpg_islands: int = 0
    category: str = ""
    sequence: str = ""


@dataclass
class ValidationResult:
    gene_name: str
    organism: str
    published_cai: float
    biocompiler_cai: float
    dnachisel_cai: float
    biocompiler_improvement: float
    dnachisel_improvement: float
    citation: str
    native_cai: float = 0.0


# ═══════════════════════════════════════════════════════════════════════
# Gene Panel Assembly
# ═══════════════════════════════════════════════════════════════════════

def _normalize_organism(org: str) -> str:
    """Normalize organism name to CODON_ADAPTIVENESS_TABLES key."""
    mapping = {
        "Escherichia coli": "Escherichia_coli",
        "Homo sapiens": "Homo_sapiens",
        "Saccharomyces cerevisiae": "Saccharomyces_cerevisiae",
        "Mus musculus": "Mus_musculus",
        "CHO K1": "CHO_K1",
    }
    return mapping.get(org, org.replace(" ", "_"))


def build_full_gene_panel() -> list[dict]:
    """Assemble the complete gene panel from all gene sets."""
    panel = []

    # E. coli Extended
    for name, entry in E_COLI_EXTENDED.items():
        panel.append({
            "name": f"{name}_Ecoli",
            "protein": entry["protein"],
            "organism": _normalize_organism(entry["organism"]),
            "category": entry.get("category", "housekeeping"),
            "gene_set": "E_COLI_EXTENDED",
        })

    # Human Therapeutic (legacy format)
    for name, entry in HUMAN_THERAPEUTIC_GENES.items():
        panel.append({
            "name": f"{name}_Human",
            "protein": entry["protein_sequence"],
            "organism": _normalize_organism(entry["organism"]),
            "category": "therapeutic",
            "gene_set": "HUMAN_THERAPEUTIC_GENES",
        })

    # Human Therapeutic (new format)
    for name, entry in HUMAN_THERAPEUTIC.items():
        if not any(p["name"] == f"{name}_Human" for p in panel):
            panel.append({
                "name": f"{name}_Human",
                "protein": entry["protein"],
                "organism": _normalize_organism(entry["organism"]),
                "category": entry.get("category", "therapeutic"),
                "gene_set": "HUMAN_THERAPEUTIC",
            })

    # Human Signaling
    for name, entry in HUMAN_SIGNALING.items():
        panel.append({
            "name": f"{name}_Human",
            "protein": entry["protein"],
            "organism": _normalize_organism(entry["organism"]),
            "category": entry.get("category", "signaling"),
            "gene_set": "HUMAN_SIGNALING",
        })

    # Yeast Industrial
    for name, entry in YEAST_INDUSTRIAL.items():
        panel.append({
            "name": f"{name}_Yeast",
            "protein": entry.get("protein", entry.get("protein_sequence", "")),
            "organism": _normalize_organism(entry.get("organism", "Saccharomyces cerevisiae")),
            "category": entry.get("category", "industrial"),
            "gene_set": "YEAST_INDUSTRIAL",
        })

    # Vaccine Antigens
    for name, entry in VACCINE_ANTIGEN_GENES.items():
        org = _normalize_organism(entry.get("organism", "Homo sapiens"))
        # Map vaccine antigens to the best available codon table
        if org not in CODON_ADAPTIVENESS_TABLES:
            if "coli" in org.lower() or "mycobacterium" in org.lower():
                org = "Escherichia_coli"
            elif "virus" in org.lower() or "coronavirus" in org.lower() or "influenza" in org.lower():
                org = "Homo_sapiens"  # Viral proteins expressed in human cells
            elif "plasmodium" in org.lower():
                org = "Homo_sapiens"
            else:
                org = "Homo_sapiens"
        panel.append({
            "name": f"{name}_Vaccine",
            "protein": entry.get("protein_sequence", entry.get("protein", "")),
            "organism": org,
            "category": "vaccine_antigen",
            "gene_set": "VACCINE_ANTIGEN_GENES",
        })

    # Stress Test
    for name, entry in STRESS_TEST_GENES.items():
        if entry.get("organism") == "synthetic":
            continue  # Skip synthetic organisms
        panel.append({
            "name": f"{name}_Stress",
            "protein": entry.get("protein_sequence", entry.get("protein", "")),
            "organism": "Escherichia_coli",  # Default for stress tests
            "category": entry.get("stress_category", "stress"),
            "gene_set": "STRESS_TEST",
        })

    # Deduplicate by name
    seen = set()
    unique_panel = []
    for entry in panel:
        if entry["name"] not in seen and entry["protein"]:
            seen.add(entry["name"])
            unique_panel.append(entry)

    return unique_panel


# ═══════════════════════════════════════════════════════════════════════
# Optimization Functions
# ═══════════════════════════════════════════════════════════════════════

def optimize_biocompiler(protein: str, organism: str) -> GeneResult:
    t0 = time.perf_counter()
    try:
        result = optimize_sequence(
            target_protein=protein,
            organism=organism,
            gc_lo=GC_LO,
            gc_hi=GC_HI,
        )
        elapsed = (time.perf_counter() - t0) * 1000
        seq = result.sequence
        cai = compute_cai_validated(seq, organism)
        gc = gc_content(seq)
        gc_profile = compute_gc_distribution(seq)
        rs_dict = count_restriction_sites(seq, list(STANDARD_ENZYME_PANEL))
        css = count_cryptic_splice_sites(seq)
        mrna = compute_mrna_stability_score(seq)
        cpg = count_cpg_islands(seq)
        return GeneResult(
            gene_name="", organism=organism, tool="biocompiler",
            protein_length=len(protein), cai=cai, gc_content=gc,
            gc_std=gc_profile.std, runtime_ms=elapsed, success=True,
            restriction_site_total=sum(rs_dict.values()),
            cryptic_splice_sites=css, mrna_stability=mrna, cpg_islands=cpg,
            sequence=seq,
        )
    except Exception as e:
        elapsed = (time.perf_counter() - t0) * 1000
        return GeneResult(
            gene_name="", organism=organism, tool="biocompiler",
            protein_length=len(protein), cai=0.0, gc_content=0.0,
            gc_std=0.0, runtime_ms=elapsed, success=False, error=str(e),
        )


def optimize_dnachisel(protein: str, organism: str) -> GeneResult:
    if not _DNACHISEL_AVAILABLE:
        return GeneResult(
            gene_name="", organism=organism, tool="dnachisel",
            protein_length=len(protein), cai=0.0, gc_content=0.0,
            gc_std=0.0, runtime_ms=0.0, success=False,
            error="DNAchisel not installed",
        )

    species = _SPECIES_MAP.get(organism, "e_coli")
    preferred = PREFERRED_CODON_TABLES.get(organism, {})
    if preferred:
        start_seq = "".join(preferred.get(aa, AA_TO_CODONS.get(aa, ["AAA"])[0]) for aa in protein)
    else:
        start_seq = "".join(AA_TO_CODONS.get(aa, ["AAA"])[0] for aa in protein)

    t0 = time.perf_counter()
    try:
        constraints = [
            EnforceTranslation(translation=protein),
            EnforceGCContent(mini=GC_LO, maxi=GC_HI),
        ]
        objectives = [CodonOptimize(species=species)]
        problem = DnaOptimizationProblem(
            sequence=start_seq,
            constraints=constraints,
            objectives=objectives,
        )
        problem.resolve_constraints()
        problem.optimize()
        elapsed = (time.perf_counter() - t0) * 1000

        seq = str(problem.sequence)
        cai = compute_cai_validated(seq, organism)
        gc = gc_content(seq)
        gc_profile = compute_gc_distribution(seq)
        rs_dict = count_restriction_sites(seq, list(STANDARD_ENZYME_PANEL))
        css = count_cryptic_splice_sites(seq)
        mrna = compute_mrna_stability_score(seq)
        cpg = count_cpg_islands(seq)
        return GeneResult(
            gene_name="", organism=organism, tool="dnachisel",
            protein_length=len(protein), cai=cai, gc_content=gc,
            gc_std=gc_profile.std, runtime_ms=elapsed, success=True,
            restriction_site_total=sum(rs_dict.values()),
            cryptic_splice_sites=css, mrna_stability=mrna, cpg_islands=cpg,
            sequence=seq,
        )
    except Exception as e:
        elapsed = (time.perf_counter() - t0) * 1000
        return GeneResult(
            gene_name="", organism=organism, tool="dnachisel",
            protein_length=len(protein), cai=0.0, gc_content=0.0,
            gc_std=0.0, runtime_ms=elapsed, success=False, error=str(e),
        )


# ═══════════════════════════════════════════════════════════════════════
# Published Data Validation
# ═══════════════════════════════════════════════════════════════════════

def validate_against_published() -> list[ValidationResult]:
    """Validate optimized sequences against published CAI values."""
    results = []
    for (gene, organism), data in VALIDATION_SEQUENCES.items():
        protein = data.get("protein_sequence", "")
        dna = data.get("dna_sequence", "")
        published_cai = data.get("expected_cai", 0.0)
        citation = data.get("citation", "")
        org_key = _normalize_organism(organism)

        # Compute native sequence CAI
        native_cai = 0.0
        if dna and org_key in CODON_ADAPTIVENESS_TABLES:
            native_cai = compute_cai_validated(dna, org_key)

        # BioCompiler optimization
        bc_cai = 0.0
        if protein and org_key in CODON_ADAPTIVENESS_TABLES:
            try:
                bc_result = optimize_sequence(
                    target_protein=protein,
                    organism=org_key,
                    gc_lo=GC_LO,
                    gc_hi=GC_HI,
                )
                bc_cai = compute_cai_validated(bc_result.sequence, org_key)
            except Exception:
                pass

        # DNAchisel optimization
        dc_cai = 0.0
        if _DNACHISEL_AVAILABLE and protein and org_key in CODON_ADAPTIVENESS_TABLES:
            try:
                species = _SPECIES_MAP.get(org_key, "e_coli")
                preferred = PREFERRED_CODON_TABLES.get(org_key, {})
                if preferred:
                    start_seq = "".join(preferred.get(aa, AA_TO_CODONS.get(aa, ["AAA"])[0]) for aa in protein)
                else:
                    start_seq = "".join(AA_TO_CODONS.get(aa, ["AAA"])[0] for aa in protein)
                constraints = [
                    EnforceTranslation(translation=protein),
                    EnforceGCContent(mini=GC_LO, maxi=GC_HI),
                ]
                objectives = [CodonOptimize(species=species)]
                problem = DnaOptimizationProblem(
                    sequence=start_seq,
                    constraints=constraints,
                    objectives=objectives,
                )
                problem.resolve_constraints()
                problem.optimize()
                dc_cai = compute_cai_validated(str(problem.sequence), org_key)
            except Exception:
                pass

        results.append(ValidationResult(
            gene_name=gene,
            organism=organism,
            published_cai=published_cai,
            biocompiler_cai=bc_cai,
            dnachisel_cai=dc_cai,
            biocompiler_improvement=bc_cai - native_cai if native_cai > 0 else bc_cai - published_cai,
            dnachisel_improvement=dc_cai - native_cai if native_cai > 0 else dc_cai - published_cai,
            citation=citation,
            native_cai=native_cai,
        ))
    return results


# ═══════════════════════════════════════════════════════════════════════
# Statistical Analysis
# ═══════════════════════════════════════════════════════════════════════

def compute_statistics(bc_values: list[float], dc_values: list[float]) -> dict:
    """Compute paired t-test, Wilcoxon signed-rank, Cohen's d."""
    import numpy as np
    from scipy import stats as sp_stats

    n = len(bc_values)
    if n < 2:
        return {"n": n, "insufficient_data": True}

    bc_arr = np.array(bc_values, dtype=np.float64)
    dc_arr = np.array(dc_values, dtype=np.float64)
    diffs = bc_arr - dc_arr

    t_stat, t_pval = sp_stats.ttest_rel(bc_arr, dc_arr)
    if np.all(diffs == 0):
        w_stat, w_pval = 0.0, 1.0
    else:
        try:
            result = sp_stats.wilcoxon(bc_arr, dc_arr, alternative="two-sided")
            w_stat, w_pval = float(result.statistic), float(result.pvalue)
        except ValueError:
            w_stat, w_pval = 0.0, 1.0

    std_diff = float(np.std(diffs, ddof=1)) if n > 1 else 0.0
    cohens_d = float(np.mean(diffs)) / std_diff if std_diff > 0 else 0.0

    bc_wins = int(np.sum(bc_arr > dc_arr + 0.001))
    dc_wins = int(np.sum(dc_arr > bc_arr + 0.001))
    ties = n - bc_wins - dc_wins

    return {
        "n": n,
        "mean_bc": float(np.mean(bc_arr)),
        "mean_dc": float(np.mean(dc_arr)),
        "mean_diff": float(np.mean(diffs)),
        "std_diff": std_diff,
        "t_statistic": float(t_stat),
        "t_pvalue": float(t_pval),
        "wilcoxon_statistic": w_stat,
        "wilcoxon_pvalue": w_pval,
        "cohens_d": cohens_d,
        "bc_wins": bc_wins,
        "dc_wins": dc_wins,
        "ties": ties,
    }


# ═══════════════════════════════════════════════════════════════════════
# Plotting
# ═══════════════════════════════════════════════════════════════════════

def plot_cai_comparison(all_results: list[GeneResult], output_dir: Path) -> None:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import matplotlib.font_manager as fm
    fm.fontManager.addfont('/usr/share/fonts/truetype/chinese/SarasaMonoSC-Regular.ttf')
    plt.rcParams['font.sans-serif'] = ['DejaVu Sans', 'Sarasa Mono SC']
    plt.rcParams['axes.unicode_minus'] = False
    import numpy as np

    organisms = sorted(set(r.organism for r in all_results if r.success))
    n_orgs = len(organisms)

    fig, axes = plt.subplots(1, max(n_orgs, 1), figsize=(7 * n_orgs, 8), sharey=True)
    if n_orgs == 1:
        axes = [axes]

    fig.suptitle("CAI Comparison: BioCompiler vs DNAchisel\n(Validated CAI — Sharp & Li 1987)", fontsize=14, fontweight="bold")

    colors_bc = "#2ecc71"
    colors_dc = "#e74c3c"

    for idx, org in enumerate(organisms):
        ax = axes[idx]
        bc_results = {r.gene_name: r for r in all_results
                      if r.organism == org and r.tool == "biocompiler" and r.success}
        dc_results = {r.gene_name: r for r in all_results
                      if r.organism == org and r.tool == "dnachisel" and r.success}
        common = sorted(set(bc_results.keys()) & set(dc_results.keys()))

        if not common:
            ax.set_title(f"{org.replace('_', ' ')}\nNo paired data")
            continue

        bc_cais = [bc_results[g].cai for g in common]
        dc_cais = [dc_results[g].cai for g in common]

        x = np.arange(len(common))
        width = 0.35
        ax.bar(x - width/2, bc_cais, width, label="BioCompiler", color=colors_bc, alpha=0.85)
        ax.bar(x + width/2, dc_cais, width, label="DNAchisel", color=colors_dc, alpha=0.85)

        ax.set_xlabel("Gene")
        ax.set_ylabel("CAI")
        org_display = org.replace("Escherichia_coli", "E. coli").replace("Homo_sapiens", "Human").replace("Saccharomyces_cerevisiae", "Yeast")
        ax.set_title(f"{org_display}\n(n={len(common)})")
        ax.set_xticks(x)
        short_names = [g[:10] for g in common]
        ax.set_xticklabels(short_names, rotation=45, ha="right", fontsize=7)
        ax.legend(loc="best", fontsize=9)
        ax.set_ylim(0, 1.05)
        ax.grid(axis="y", alpha=0.3)

    plt.tight_layout()
    fig.savefig(output_dir / "cai_comparison.png", dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"  CAI comparison plot saved")


def plot_speed_comparison(all_results: list[GeneResult], output_dir: Path) -> None:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import matplotlib.font_manager as fm
    fm.fontManager.addfont('/usr/share/fonts/truetype/chinese/SarasaMonoSC-Regular.ttf')
    plt.rcParams['font.sans-serif'] = ['DejaVu Sans', 'Sarasa Mono SC']
    plt.rcParams['axes.unicode_minus'] = False
    import numpy as np

    fig, ax = plt.subplots(figsize=(10, 7))
    fig.suptitle("Runtime Comparison: BioCompiler vs DNAchisel", fontsize=14, fontweight="bold")

    bc_lengths, bc_times = [], []
    dc_lengths, dc_times = [], []

    organisms = sorted(set(r.organism for r in all_results if r.success))
    for org in organisms:
        bc_dict = {r.gene_name: r for r in all_results
                   if r.organism == org and r.tool == "biocompiler" and r.success}
        dc_dict = {r.gene_name: r for r in all_results
                   if r.organism == org and r.tool == "dnachisel" and r.success}
        common = set(bc_dict.keys()) & set(dc_dict.keys())
        for g in common:
            bc_lengths.append(bc_dict[g].protein_length)
            bc_times.append(bc_dict[g].runtime_ms)
            dc_lengths.append(dc_dict[g].protein_length)
            dc_times.append(dc_dict[g].runtime_ms)

    ax.scatter(bc_lengths, bc_times, c="#2ecc71", s=70, alpha=0.7, label="BioCompiler", zorder=3, edgecolors="darkgreen", linewidth=0.5)
    ax.scatter(dc_lengths, dc_times, c="#e74c3c", s=70, alpha=0.7, marker="^", label="DNAchisel", zorder=3, edgecolors="darkred", linewidth=0.5)

    ax.set_xlabel("Protein Length (amino acids)")
    ax.set_ylabel("Runtime (ms)")
    ax.legend(loc="best")
    ax.grid(alpha=0.3)
    ax.set_yscale("log")

    plt.tight_layout()
    fig.savefig(output_dir / "speed_comparison.png", dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"  Speed comparison plot saved")


def plot_published_validation(validation_results: list[ValidationResult], output_dir: Path) -> None:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import matplotlib.font_manager as fm
    fm.fontManager.addfont('/usr/share/fonts/truetype/chinese/SarasaMonoSC-Regular.ttf')
    plt.rcParams['font.sans-serif'] = ['DejaVu Sans', 'Sarasa Mono SC']
    plt.rcParams['axes.unicode_minus'] = False
    import numpy as np

    fig, ax = plt.subplots(figsize=(14, 7))
    fig.suptitle("Published CAI Validation: Optimized vs Native vs Published", fontsize=14, fontweight="bold")

    labels = [f"{v.gene_name}\n({v.organism.replace('Escherichia_coli','E.coli').replace('Homo_sapiens','Human').replace('Saccharomyces_cerevisiae','Yeast')})"
              for v in validation_results]
    published = [v.published_cai for v in validation_results]
    native = [v.native_cai for v in validation_results]
    bc_cai = [v.biocompiler_cai for v in validation_results]
    dc_cai = [v.dnachisel_cai for v in validation_results]

    x = np.arange(len(labels))
    width = 0.2

    ax.bar(x - 1.5*width, published, width, label="Published (Literature)", color="#3498db", alpha=0.8)
    ax.bar(x - 0.5*width, native, width, label="Native (RefSeq)", color="#95a5a6", alpha=0.8)
    ax.bar(x + 0.5*width, bc_cai, width, label="BioCompiler Optimized", color="#2ecc71", alpha=0.85)
    ax.bar(x + 1.5*width, dc_cai, width, label="DNAchisel Optimized", color="#e74c3c", alpha=0.85)

    ax.set_xlabel("Gene / Organism")
    ax.set_ylabel("CAI")
    ax.set_xticks(x)
    ax.set_xticklabels(labels, rotation=45, ha="right", fontsize=7)
    ax.legend(loc="best", fontsize=9)
    ax.set_ylim(0, 1.05)
    ax.grid(axis="y", alpha=0.3)

    plt.tight_layout()
    fig.savefig(output_dir / "published_validation.png", dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"  Published validation plot saved")


def plot_cai_by_category(all_results: list[GeneResult], output_dir: Path) -> None:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import matplotlib.font_manager as fm
    fm.fontManager.addfont('/usr/share/fonts/truetype/chinese/SarasaMonoSC-Regular.ttf')
    plt.rcParams['font.sans-serif'] = ['DejaVu Sans', 'Sarasa Mono SC']
    plt.rcParams['axes.unicode_minus'] = False
    import numpy as np

    categories = sorted(set(r.category for r in all_results if r.success and r.category))

    fig, ax = plt.subplots(figsize=(12, 6))
    fig.suptitle("CAI by Gene Category: BioCompiler vs DNAchisel", fontsize=14, fontweight="bold")

    bc_means, dc_means, bc_stds, dc_stds, cat_labels = [], [], [], [], []
    for cat in categories:
        bc_cais = [r.cai for r in all_results if r.category == cat and r.tool == "biocompiler" and r.success]
        dc_cais = [r.cai for r in all_results if r.category == cat and r.tool == "dnachisel" and r.success]
        if bc_cais and dc_cais:
            bc_means.append(np.mean(bc_cais))
            dc_means.append(np.mean(dc_cais))
            bc_stds.append(np.std(bc_cais))
            dc_stds.append(np.std(dc_cais))
            cat_labels.append(cat)

    x = np.arange(len(cat_labels))
    width = 0.35
    ax.bar(x - width/2, bc_means, width, yerr=bc_stds, label="BioCompiler", color="#2ecc71", alpha=0.85, capsize=3)
    ax.bar(x + width/2, dc_means, width, yerr=dc_stds, label="DNAchisel", color="#e74c3c", alpha=0.85, capsize=3)

    ax.set_xlabel("Gene Category")
    ax.set_ylabel("Mean CAI")
    ax.set_xticks(x)
    ax.set_xticklabels(cat_labels, rotation=30, ha="right", fontsize=9)
    ax.legend(loc="best")
    ax.set_ylim(0, 1.05)
    ax.grid(axis="y", alpha=0.3)

    plt.tight_layout()
    fig.savefig(output_dir / "cai_by_category.png", dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"  CAI by category plot saved")


# ═══════════════════════════════════════════════════════════════════════
# Main Benchmark Runner
# ═══════════════════════════════════════════════════════════════════════

def run_mega_benchmark(output_dir: Path) -> None:
    output_dir.mkdir(parents=True, exist_ok=True)

    print("=" * 90)
    print("  MEGA BENCHMARK: BioCompiler vs DNAchisel — Full Gene Panel")
    print("  CAI computed with compute_cai_validated for BOTH tools (Sharp & Li 1987)")
    print("=" * 90)
    print()

    if _DNACHISEL_AVAILABLE:
        print("  DNAchisel: AVAILABLE (CodonOptimize objective + validated CAI)")
    else:
        print("  DNAchisel: NOT INSTALLED — BioCompiler-only results")
    print()

    # ── Build gene panel ──
    gene_panel = build_full_gene_panel()
    print(f"  Gene panel: {len(gene_panel)} genes across {len(set(g['organism'] for g in gene_panel))} organisms")
    for org in sorted(set(g['organism'] for g in gene_panel)):
        count = sum(1 for g in gene_panel if g['organism'] == org)
        print(f"    {org}: {count} genes")
    print()

    # ── Run head-to-head benchmark ──
    all_results: list[GeneResult] = []

    print(f"  {'Gene':<25} {'Organism':<22} {'Len':>4} │ "
          f"{'BC CAI':>8} {'DC CAI':>8} {'ΔCAI':>8} │ "
          f"{'BC ms':>8} {'DC ms':>8} {'Spd Δ':>8}")
    print("  " + "-" * 110)

    for entry in gene_panel:
        name = entry["name"]
        protein = entry["protein"]
        organism = entry["organism"]
        category = entry.get("category", "")
        aa_len = len(protein)

        # Skip if organism not supported
        if organism not in CODON_ADAPTIVENESS_TABLES:
            continue

        # BioCompiler
        bc = optimize_biocompiler(protein, organism)
        bc.gene_name = name
        bc.category = category
        all_results.append(bc)

        # DNAchisel
        dc = optimize_dnachisel(protein, organism)
        dc.gene_name = name
        dc.category = category
        all_results.append(dc)

        # Print inline
        bc_cai = bc.cai if bc.success else 0.0
        dc_cai = dc.cai if dc.success else 0.0
        delta_cai = bc_cai - dc_cai
        bc_ms = bc.runtime_ms
        dc_ms = dc.runtime_ms if dc.success else 0.0
        speed_delta = f"{dc_ms/bc_ms:.1f}x" if bc_ms > 0 and dc.success else "N/A"

        err_flag = ""
        if not bc.success:
            err_flag += " [BC ERR]"
        if not dc.success:
            err_flag += " [DC N/A]"

        org_display = organism.replace("Escherichia_coli", "E.coli").replace("Homo_sapiens", "Human").replace("Saccharomyces_cerevisiae", "Yeast")
        print(f"  {name:<25} {org_display:<22} {aa_len:>4} │ "
              f"{bc_cai:>8.4f} {dc_cai:>8.4f} {delta_cai:>+8.4f} │ "
              f"{bc_ms:>8.1f} {dc_ms:>8.1f} {speed_delta:>8}{err_flag}")

    print("  " + "-" * 110)
    print()

    # ── Published data validation ──
    print("  VALIDATION AGAINST PUBLISHED DATA")
    print("  " + "-" * 70)
    validation_results = validate_against_published()

    for v in validation_results:
        org_display = v.organism.replace("Escherichia_coli", "E.coli").replace("Homo_sapiens", "Human").replace("Saccharomyces_cerevisiae", "Yeast")
        bc_imp = f"+{v.biocompiler_improvement:.3f}" if v.biocompiler_improvement > 0 else f"{v.biocompiler_improvement:.3f}"
        dc_imp = f"+{v.dnachisel_improvement:.3f}" if v.dnachisel_improvement > 0 else f"{v.dnachisel_improvement:.3f}"
        print(f"  {v.gene_name:<12} {org_display:<10} │ Pub={v.published_cai:.2f}  Nat={v.native_cai:.3f}  "
              f"BC={v.biocompiler_cai:.4f}({bc_imp})  DC={v.dnachisel_cai:.4f}({dc_imp})")

    print("  " + "-" * 70)
    print()

    # ── Statistical analysis per organism ──
    print("  STATISTICAL ANALYSIS BY ORGANISM")
    print("  " + "-" * 90)

    organisms_with_data = sorted(set(r.organism for r in all_results if r.success))
    all_stats = {}

    for org in organisms_with_data:
        bc_results = {r.gene_name: r for r in all_results
                      if r.organism == org and r.tool == "biocompiler" and r.success}
        dc_results = {r.gene_name: r for r in all_results
                      if r.organism == org and r.tool == "dnachisel" and r.success}
        common = sorted(set(bc_results.keys()) & set(dc_results.keys()))

        if len(common) < 2:
            print(f"  {org}: Insufficient paired data (n={len(common)})")
            continue

        # CAI stats
        bc_cais = [bc_results[g].cai for g in common]
        dc_cais = [dc_results[g].cai for g in common]
        cai_stats = compute_statistics(bc_cais, dc_cais)

        # Speed stats
        bc_times = [bc_results[g].runtime_ms for g in common]
        dc_times = [dc_results[g].runtime_ms for g in common]
        speed_stats = compute_statistics(bc_times, dc_times)

        all_stats[org] = {"cai": cai_stats, "speed": speed_stats}

        org_display = org.replace("Escherichia_coli", "E. coli").replace("Homo_sapiens", "Human").replace("Saccharomyces_cerevisiae", "Yeast")
        print(f"\n  {org_display} (n={len(common)} genes):")
        print(f"    CAI:    BC={cai_stats['mean_bc']:.4f}  DC={cai_stats['mean_dc']:.4f}  "
              f"Δ={cai_stats['mean_diff']:+.4f}  t p={cai_stats['t_pvalue']:.4f}  "
              f"W p={cai_stats['wilcoxon_pvalue']:.4f}  d={cai_stats['cohens_d']:.2f}  "
              f"Wins: BC={cai_stats['bc_wins']} DC={cai_stats['dc_wins']} Ties={cai_stats['ties']}")
        print(f"    Speed:  BC={speed_stats['mean_bc']:.1f}ms  DC={speed_stats['mean_dc']:.1f}ms  "
              f"Ratio={speed_stats['mean_bc']/speed_stats['mean_dc']:.2f}x  "
              f"t p={speed_stats['t_pvalue']:.4f}  d={speed_stats['cohens_d']:.2f}")

    print("  " + "-" * 90)
    print()

    # ── Rank correlation with expression level ──
    print("  EXPRESSION LEVEL RANK CORRELATION (E. coli)")
    print("  " + "-" * 60)
    # E. coli genes have expected_cai_range that correlates with expression
    ecoli_gene_names = {name for name, _ in [(k, None) for k in E_COLI_EXTENDED.keys()]}
    bc_ecoli = {r.gene_name: r for r in all_results
                if r.organism == "Escherichia_coli" and r.tool == "biocompiler" and r.success
                and r.gene_name.replace("_Ecoli", "") in ecoli_gene_names}
    dc_ecoli = {r.gene_name: r for r in all_results
                if r.organism == "Escherichia_coli" and r.tool == "dnachisel" and r.success
                and r.gene_name.replace("_Ecoli", "") in ecoli_gene_names}

    # Expected CAI ranges from gene_sets (higher = more highly expressed)
    expected_ranks = {}
    for name, entry in E_COLI_EXTENDED.items():
        expected_ranks[f"{name}_Ecoli"] = entry.get("expected_cai_range", (0.5, 0.7))[1]  # Use upper bound

    common_ecoli = sorted(set(bc_ecoli.keys()) & set(dc_ecoli.keys()) & set(expected_ranks.keys()))
    if len(common_ecoli) >= 3:
        from scipy import stats as sp_stats
        import numpy as np

        expected_vals = [expected_ranks[g] for g in common_ecoli]
        bc_cais_rank = [bc_ecoli[g].cai for g in common_ecoli]
        dc_cais_rank = [dc_ecoli[g].cai for g in common_ecoli]

        bc_spearman = sp_stats.spearmanr(expected_vals, bc_cais_rank)
        dc_spearman = sp_stats.spearmanr(expected_vals, dc_cais_rank)

        print(f"    BioCompiler Spearman ρ = {bc_spearman.correlation:.3f} (p={bc_spearman.pvalue:.4f})")
        print(f"    DNAchisel   Spearman ρ = {dc_spearman.correlation:.3f} (p={dc_spearman.pvalue:.4f})")
        print(f"    (Higher ρ = better rank correlation with expression level)")
    print()

    # ── Save JSON results ──
    json_data = {
        "benchmark_metadata": {
            "timestamp": time.strftime("%Y-%m-%d %H:%M:%S"),
            "dnachisel_available": _DNACHISEL_AVAILABLE,
            "gc_range": [GC_LO, GC_HI],
            "total_genes": len(gene_panel),
        },
        "head_to_head_results": [],
        "validation_results": [],
        "statistics": {org: stats for org, stats in all_stats.items()},
    }

    for r in all_results:
        d = asdict(r)
        d["cai"] = round(d["cai"], 6)
        d["gc_content"] = round(d["gc_content"], 6)
        d["runtime_ms"] = round(d["runtime_ms"], 3)
        del d["sequence"]  # Too large for JSON
        json_data["head_to_head_results"].append(d)

    for v in validation_results:
        json_data["validation_results"].append(asdict(v))

    json_path = output_dir / "mega_benchmark_results.json"
    with open(json_path, "w") as f:
        json.dump(json_data, f, indent=2)
    print(f"  JSON results saved to {json_path}")

    # ── Save CSV ──
    csv_path = output_dir / "mega_benchmark_results.csv"
    with open(csv_path, "w") as f:
        headers = ["gene_name", "organism", "tool", "protein_length", "cai", "gc_content",
                   "gc_std", "runtime_ms", "success", "restriction_site_total",
                   "cryptic_splice_sites", "mrna_stability", "cpg_islands", "category"]
        f.write(",".join(headers) + "\n")
        for r in all_results:
            row = [
                r.gene_name, r.organism, r.tool, r.protein_length,
                f"{r.cai:.6f}", f"{r.gc_content:.6f}", f"{r.gc_std:.6f}",
                f"{r.runtime_ms:.3f}", r.success, r.restriction_site_total,
                r.cryptic_splice_sites, f"{r.mrna_stability:.4f}", r.cpg_islands,
                r.category,
            ]
            f.write(",".join(str(v) for v in row) + "\n")
    print(f"  CSV results saved to {csv_path}")

    # ── Generate plots ──
    try:
        plot_cai_comparison(all_results, output_dir)
    except Exception as e:
        print(f"  Warning: CAI plot failed: {e}")

    try:
        plot_speed_comparison(all_results, output_dir)
    except Exception as e:
        print(f"  Warning: Speed plot failed: {e}")

    try:
        plot_published_validation(validation_results, output_dir)
    except Exception as e:
        print(f"  Warning: Published validation plot failed: {e}")

    try:
        plot_cai_by_category(all_results, output_dir)
    except Exception as e:
        print(f"  Warning: Category plot failed: {e}")

    # ── Overall Summary ──
    print()
    print("=" * 90)
    print("  OVERALL SUMMARY")
    print("=" * 90)

    total_bc_success = sum(1 for r in all_results if r.tool == "biocompiler" and r.success)
    total_dc_success = sum(1 for r in all_results if r.tool == "dnachisel" and r.success)

    bc_cais_all = [r.cai for r in all_results if r.tool == "biocompiler" and r.success]
    dc_cais_all = [r.cai for r in all_results if r.tool == "dnachisel" and r.success]
    bc_times_all = [r.runtime_ms for r in all_results if r.tool == "biocompiler" and r.success]
    dc_times_all = [r.runtime_ms for r in all_results if r.tool == "dnachisel" and r.success]

    if bc_cais_all:
        print(f"  BioCompiler: {total_bc_success} genes optimized")
        print(f"    Mean CAI: {sum(bc_cais_all)/len(bc_cais_all):.4f}")
        print(f"    Mean Time: {sum(bc_times_all)/len(bc_times_all):.1f} ms")

    if dc_cais_all:
        print(f"  DNAchisel: {total_dc_success} genes optimized")
        print(f"    Mean CAI: {sum(dc_cais_all)/len(dc_cais_all):.4f}")
        print(f"    Mean Time: {sum(dc_times_all)/len(dc_times_all):.1f} ms")

    if bc_cais_all and dc_cais_all:
        overall_cai_gap = sum(bc_cais_all)/len(bc_cais_all) - sum(dc_cais_all)/len(dc_cais_all)
        overall_speed_ratio = (sum(bc_times_all)/len(bc_times_all)) / (sum(dc_times_all)/len(dc_times_all)) if dc_times_all else 0
        print(f"\n  CAI Gap (BC - DC): {overall_cai_gap:+.4f}")
        print(f"  Speed Ratio (BC/DC): {overall_speed_ratio:.2f}x")
        if overall_cai_gap >= 0:
            print(f"\n  >>> BioCompiler MATCHES or EXCEEDS DNAchisel on CAI <<<")
        if overall_speed_ratio < 1.0:
            print(f"  >>> BioCompiler is FASTER than DNAchisel <<<")
        elif overall_speed_ratio < 2.0:
            print(f"  >>> BioCompiler is within 2x of DNAchisel on speed <<<")

    # Validation summary
    if validation_results:
        bc_above_pub = sum(1 for v in validation_results if v.biocompiler_cai > v.published_cai)
        dc_above_pub = sum(1 for v in validation_results if v.dnachisel_cai > v.published_cai)
        n_val = len(validation_results)
        print(f"\n  Published CAI Validation: BC above published {bc_above_pub}/{n_val}, DC above published {dc_above_pub}/{n_val}")

    print("=" * 90)


if __name__ == "__main__":
    output_dir = Path("/home/z/my-project/download/benchmark_results")
    run_mega_benchmark(output_dir)
