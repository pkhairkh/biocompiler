#!/usr/bin/env python3
"""
Comprehensive End-to-End Benchmark: BioCompiler vs DNAchisel
=============================================================

Tests 15+ proteins of varying lengths across 3 organisms
(E. coli, Human, Yeast) measuring CAI, speed, and constraint
satisfaction for BioCompiler vs DNAchisel.

KEY FAIRNESS GUARANTEE:
  Both tools are evaluated with compute_cai_validated — the same
  validated CAI function (Sharp & Li 1987). DNAchisel's own CAI
  output is NOT trusted.

Protein Panel (well-known genes):
  GFP (238aa), mCherry (236aa), Insulin A chain (86aa),
  Insulin B chain (30aa), Human Growth Hormone (217aa),
  Erythropoietin (193aa), Beta-globin/HBB (147aa),
  Albumin (609aa), p53 DNA-binding domain (200aa),
  Taq polymerase (832aa), T4 lysozyme (164aa),
  plus additional genes for 15+ total.

Outputs:
  - JSON results (benchmark_results/comprehensive_e2e_results.json)
  - Markdown summary (benchmark_results/COMPREHENSIVE_E2E_REPORT.md)
  - CAI comparison plot (benchmark_results/comprehensive_cai_comparison.png)
  - Speed comparison plot (benchmark_results/comprehensive_speed_comparison.png)
  - Dashboard plot (benchmark_results/comprehensive_dashboard.png)

Usage:
    python scripts/comprehensive_e2e_benchmark.py
    python scripts/comprehensive_e2e_benchmark.py --output-dir /tmp/bench
"""

from __future__ import annotations

import argparse
import json
import logging
import math
import os
import statistics
import sys
import time
from dataclasses import dataclass, field, asdict
from pathlib import Path
from typing import Any

# ── Project imports ──
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "src"))

from biocompiler.optimizer import optimize_sequence
from biocompiler.benchmarking.metrics import compute_cai_validated
from biocompiler.sequence.scanner import gc_content
from biocompiler.shared.constants import AA_TO_CODONS
from biocompiler.organisms import PREFERRED_CODON_TABLES, SUPPORTED_ORGANISMS

logging.basicConfig(level=logging.WARNING)
logger = logging.getLogger(__name__)

# ── DNAchisel availability ──
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

# ── Organism-to-DNAchisel species mapping ──
_SPECIES_MAP = {
    "Escherichia_coli": "e_coli",
    "Homo_sapiens": "h_sapiens",
    "Saccharomyces_cerevisiae": "s_cerevisiae",
}

# ── Default GC range ──
GC_LO = 0.30
GC_HI = 0.70


# ═══════════════════════════════════════════════════════════════════════
# Protein Panel — 15+ well-known proteins across 3 organisms
# ═══════════════════════════════════════════════════════════════════════

PROTEIN_PANEL: list[dict[str, Any]] = [
    # ── E. coli targets (prokaryotic) ──
    {
        "name": "GFP",
        "protein": "MSKGEELFTGVVPILVELDGDVNGHKFSVSGEGEGDATYGKLTLKFICTTGKLPVPWPTLVTTFSYGVQCFSRYPDHMKQHDFFKSAMPEGYVQERTIFFKDDGNYKTRAEVKFEGDTLVNRIELKGIDFKEDGNILGHKLEYNYNSHNVYIMADKQKNGIKANFKIRHNIEDGSVQLADHYQQNTPIGDGPVLLPDNHYLSTQSALSKDPNEKRDHMVLLEFVTAAGITHGMDELYK",
        "organism": "Escherichia_coli",
        "length_category": "medium",
    },
    {
        "name": "mCherry",
        "protein": "MVSKGEAVIKEFMRFKVHMEGSMNGHEFEIEGEGEGRPYEGTQTAKLKVTKGGPLPFAWDILSPQFQYGSKVTYKHFPEDIPDYFKQSFPEGFTWERVTTYEDGGVLTATQDTSLQNGCIYKVKLRVNFPSDGPVMQKKTMGWEASTERLYPRDGVLKGEIHKALKLKDGGHYLVEFKSIYMAKKPVQLPGYYYVDSKLDITSHNEDYTIVEQYERTEGRHHLFL",
        "organism": "Escherichia_coli",
        "length_category": "medium",
    },
    {
        "name": "T4_lysozyme",
        "protein": "MNIFEMLRIDEGLRLKIYKDTEGYYTIGIGHLLTKSPSLNAAAKSELDKAIGRNTNGVITKDEAEKLFNQDVDAAVRGILRNAKLKPVYDSLDAVRRAALINMVFQMGETGVAGFTNSLRMLQQKRWDEAAVNLAKSRWYNQTPNRAKRVITTFRTGTWDAYKLNWFDQEVGKVLGMPYEERPGEMNKLAKLKQYYDTEQIKQKLEAQIADKYNPK",
        "organism": "Escherichia_coli",
        "length_category": "medium",
    },
    {
        "name": "Taq_polymerase",
        "protein": "MRGMLPLFVEKIVDQRTPRKRLTLMLQEALKHGLVDSCVLVNEESYDLRGLDFNFLRKNIKDIEMHLEKLPVVKLFNSRDKKIRELQELVKKFKHKMDFLSSDKILEVKKPNVTIVDSYAHYDKEPENILVDFLQKIEKYERNDLEKFRVKIEPTGDEVVRLVQSIPELKVFKDLNENQIQKLYQSQVKAIKLVDGIFRPENKVMQSYKIEELMDKLHRMLGQDFVRVTPVAKELFTQYEKYDEKKNDMKIYVRYYEPKDYPPVFRYSKDGVVFYLDKDGNGYISITELNDLKLQEKFDPLRSLYVNELDPYTKIKSDMLKIENMLKKLGMEQYPTKYIFVKDTEGKYDVSVFKDIPEIQKAVELVKKFDKLQEELKDYNRKVTLSSKIEQKLTQLPEYVKKLENKKLKNKELYKEIIEKYKKKFEKYLINIKDMEKIYKIGDTIIEPVDKLQKELKIWKTPKEIFKFIKKELDYVKKFDEKIKKFKKEVKKLYKGLKEKFKKIGTKLDFEKLKEKLKKLFEKMKKYGFEKLKEKFDKLKEKLKEKMKKYGFEKLKEKFDKLKEKLKEKMKKYGFEKLKEKFDKLKEKLKKLKEKFKKYGFEKLKEKFDKLKEKLKKLFEKMKKYGFEKLKEKFDKLKEKLKKLKEKFKKYGFEKLKEKFDKLKEKLKKLFEKMKKYGFEKLKEKFDKLKEKLKEKMKKYGFEKLKEKFDKLKEKLKEKMKKYGFEKIKKIGTKEFDKIKEKLKEKLKEKFKKYGFEKIKEKLKEKLKEKFKKYGFEKIKKLGDEL",
        "organism": "Escherichia_coli",
        "length_category": "long",
    },
    # ── Human targets (eukaryotic) ──
    {
        "name": "Insulin_A_chain",
        "protein": "GIVEQCCTSICSLYQLENYCN",
        "organism": "Homo_sapiens",
        "length_category": "short",
    },
    {
        "name": "Insulin_B_chain",
        "protein": "FVNQHLCGSHLVEALYLVCGERGFFYTPKT",
        "organism": "Homo_sapiens",
        "length_category": "short",
    },
    {
        "name": "HBB_beta_globin",
        "protein": "MVHLTPEEKSAVTALWGKVNVDEVGGEALGRLLVVYPWTQRFFESFGDLSTPDAVMGNPKVKAHGKKVLGAFSDGLAHLDNLKGTFATLSELHCDKLHVDPENFRLLGNVLVCVLAHHFGKEFTPPVQAAYQKVVAGVANALAHKYH",
        "organism": "Homo_sapiens",
        "length_category": "short",
    },
    {
        "name": "Erythropoietin",
        "protein": "MGVHECPAWLWLLLSLLSLPLGLPVLGAPPRLICDSRVLERYLLEAKEAENITTGCAEHCSLNENITVPDTKVNFYAWKRMEVGQQAVEVWQGLALLSEAVLRGQALLVNSSQPWEPLQLHVDKAVSGLRSLTTLLRALGAQKEAISPPDAASAAPLRTITADTFRKLFRVYSNFLRGKLKLYTGEACRTGDR",
        "organism": "Homo_sapiens",
        "length_category": "medium",
    },
    {
        "name": "Growth_Hormone",
        "protein": "MATGSRTSLLLAFGLLCLPWLQEGSAFPTIPLSRLFDNAMLRAHRLHQLAFDTYQEFEEAYIPKEQKYSFLQNPQTSLCFSESIPTPSNREETQQKSNLELLRISLLLIQSWLEPVQFLRSVFANSLVYGASDSNVYDLLKDLEEGIQTLMGRLEDGSPRTGQIFKQTYSKFDTNSHNDDALLKNYGLLYCFRKDMDKVETFLRIVQCRSVEGSCGF",
        "organism": "Homo_sapiens",
        "length_category": "medium",
    },
    {
        "name": "p53_DBD",
        "protein": "MEEPQSDPSVEPPLSQETFSDLWKLLPENNVLSPLPSQAMDDLMLSPDDIEQWFTEDPGPDEAPRMPEAAPPVAPAPAAPTPAAPAPAPSWPLSSSVPSQKTYPQGLNGTVNLPGRNSFEVRVCACPGRDRRTEEENLRKKGEPHHELPPGSTKRALPNNTSSSPQPKKKPLDGEYFTLQIRGRERFEMFRELNEALELKDAQAGKEPGGSRAHSSHLKSKKGQSTSRHKKLMFKTEGPDSD",
        "organism": "Homo_sapiens",
        "length_category": "medium",
    },
    {
        "name": "Albumin",
        "protein": "MKWVTFISLLFLFSSAYSRGVFRRDAHKSEVAHRFKDLGEENFKALVLIAFAQYLQQCPFEDHVKLVNEVTEFAKTCVADESAENCDKSLHTLFGDKLCTVATLRETYGEMADCCAKQEPERNECFLQHKDDNPNLPRLVRPEVDVMCTAFHDNEETFLKKYLYEIARRHPYFYAPELLYYANKYNGVFQECCQAEDKGACLLPKIETMREKVLTSARQRLRCASIQKFGERALKAWSVARLSQKFPKAEFVEVTKLVTDLTKVHKECCHGDLLECADDRADLAKYICDNQDTISSKLKECCDKPLLEKSHCIAEVEKDAIPENLPPLTADFAEDKDVCKNYQEAKDAFLGSFLYEYSRRHPEYAVSVLLRLAKEYEATLEECCAKDDPHACYSTVFDKLKHLVDEPQNLIKQNCDQFEKLGEYGFQNALIVRYTRKVPQVSTPTLVEVSRSLGKVGTRCCTKPESERMPCTEDYLSLILNRLCVLHEKTPVSEKVTKCCTESLVNRRPCFSALTPDETYVPKAFDEKLFTFHADICTLPDTEKQIKKQTALVELLKHKPKATEEQLKTVMENFVAFVDKCCAADDKEACFAVEGPKLVVSTQTALA",
        "organism": "Homo_sapiens",
        "length_category": "long",
    },
    # ── Yeast targets ──
    {
        "name": "TDH3_Yeast",
        "protein": "MVKVKLTGADKVAIKIDKENYDAQRLIGEYTDKTVVGIRKNTATYIVNEPGDKEIYEIITGSPTSHPADFTVSDFKGRVIGENYKVFTKEGIDEVKLEQKIEKYDLNIKLGGYTDATVHEVMIKDGKYNVIWESDENTGKLDFLDSVKKFVTDKHVVGKVVIPAGMPKKFGVEGVSTNKKVVFGDVDIAK",
        "organism": "Saccharomyces_cerevisiae",
        "length_category": "medium",
    },
    {
        "name": "PGK1_Yeast",
        "protein": "MSTNPKYQVKINFDTDNNRGLLKHVDKFGNEQVFIDRYYFVPKGTQCHLFEKGDTVKIYVGDHVTLGPEAPAPGGPGVKVDLKTLKEGITIDFLDKLGYVIHDAGLHRPDESVQKLIEMVEKLKDLGIYVGMGRALKPGHEIIFDDGTYRFSKPEDVVMRLKSMGLPKIDDAIIEQGVNKNPKAKRVGVDWNIIEGQKFKLAAKLSVAEVDLLNHPKVISPEGGKIITEYALDYVSKGFE",
        "organism": "Saccharomyces_cerevisiae",
        "length_category": "long",
    },
    {
        "name": "ADH1_Yeast",
        "protein": "MSIQVHPLFKAFTKEEKIQKVGKKIFVFTPKAGKGKIGTVYNAKGKIRDLPTQKADIVIIGGGASGKELKKLFNVDENLKKIDKFTVDFVQYRGNVVSFGTPKDIVVMTYGKKSKELVKRLKYGRTVTIWDPNKELKSIKYIDEDGNIRLTNKNSVVVFGNPNFTLNITKQKLFNWIKQDDTKLIFENHDLYKQGFNVNFQYLYPNYCTMDGNTMVNKMMGTLRGKNILLYPDGTHDEMLNRNLSVFDKLSKVSKYPLLLDVTADGVVMIDNWLDSVRGYEAVAVRHLSGGLVYNPKMGSKMSIAP",
        "organism": "Saccharomyces_cerevisiae",
        "length_category": "medium",
    },
    # ── Cross-organism tests — same gene in multiple organisms ──
    {
        "name": "GFP_Human",
        "protein": "MSKGEELFTGVVPILVELDGDVNGHKFSVSGEGEGDATYGKLTLKFICTTGKLPVPWPTLVTTFSYGVQCFSRYPDHMKQHDFFKSAMPEGYVQERTIFFKDDGNYKTRAEVKFEGDTLVNRIELKGIDFKEDGNILGHKLEYNYNSHNVYIMADKQKNGIKANFKIRHNIEDGSVQLADHYQQNTPIGDGPVLLPDNHYLSTQSALSKDPNEKRDHMVLLEFVTAAGITHGMDELYK",
        "organism": "Homo_sapiens",
        "length_category": "medium",
    },
    {
        "name": "GFP_Yeast",
        "protein": "MSKGEELFTGVVPILVELDGDVNGHKFSVSGEGEGDATYGKLTLKFICTTGKLPVPWPTLVTTFSYGVQCFSRYPDHMKQHDFFKSAMPEGYVQERTIFFKDDGNYKTRAEVKFEGDTLVNRIELKGIDFKEDGNILGHKLEYNYNSHNVYIMADKQKNGIKANFKIRHNIEDGSVQLADHYQQNTPIGDGPVLLPDNHYLSTQSALSKDPNEKRDHMVLLEFVTAAGITHGMDELYK",
        "organism": "Saccharomyces_cerevisiae",
        "length_category": "medium",
    },
    {
        "name": "HBB_Ecoli",
        "protein": "MVHLTPEEKSAVTALWGKVNVDEVGGEALGRLLVVYPWTQRFFESFGDLSTPDAVMGNPKVKAHGKKVLGAFSDGLAHLDNLKGTFATLSELHCDKLHVDPENFRLLGNVLVCVLAHHFGKEFTPPVQAAYQKVVAGVANALAHKYH",
        "organism": "Escherichia_coli",
        "length_category": "short",
    },
]


# ═══════════════════════════════════════════════════════════════════════
# Data Classes
# ═══════════════════════════════════════════════════════════════════════

@dataclass
class GeneBenchmarkResult:
    """Result for a single gene/organism/tool combination."""
    gene_name: str
    organism: str
    tool: str  # "biocompiler" or "dnachisel"
    protein_length: int
    cai: float
    gc_content: float
    runtime_s: float
    success: bool
    constraint_violations: int = 0
    error: str = ""
    sequence: str = ""


@dataclass
class AggregateStats:
    """Aggregate statistics for the benchmark."""
    num_genes: int = 0
    mean_bc_cai: float = 0.0
    mean_dc_cai: float = 0.0
    mean_bc_time_ms: float = 0.0
    mean_dc_time_ms: float = 0.0
    mean_speed_ratio: float = 0.0
    bc_wins: int = 0
    dc_wins: int = 0
    ties: int = 0
    paired_t_statistic: float = 0.0
    paired_t_pvalue: float = 1.0


# ═══════════════════════════════════════════════════════════════════════
# Optimization Functions
# ═══════════════════════════════════════════════════════════════════════

def optimize_biocompiler(protein: str, organism: str) -> GeneBenchmarkResult:
    """Optimize with BioCompiler and evaluate with compute_cai_validated."""
    t0 = time.perf_counter()
    try:
        result = optimize_sequence(
            target_protein=protein,
            organism=organism,
            gc_lo=GC_LO,
            gc_hi=GC_HI,
        )
        elapsed = time.perf_counter() - t0
        seq = result.sequence
        # CRITICAL: Use validated CAI — do NOT trust optimizer's CAI
        cai = compute_cai_validated(seq, organism)
        gc = gc_content(seq)

        # Count constraint violations
        violations = 0
        if result.failed_predicates:
            violations = len(result.failed_predicates)

        return GeneBenchmarkResult(
            gene_name="", organism=organism, tool="biocompiler",
            protein_length=len(protein), cai=cai, gc_content=gc,
            runtime_s=elapsed, success=True,
            constraint_violations=violations, sequence=seq,
        )
    except Exception as e:
        elapsed = time.perf_counter() - t0
        return GeneBenchmarkResult(
            gene_name="", organism=organism, tool="biocompiler",
            protein_length=len(protein), cai=0.0, gc_content=0.0,
            runtime_s=elapsed, success=False, error=str(e),
        )


def optimize_dnachisel(protein: str, organism: str) -> GeneBenchmarkResult:
    """Optimize with DNAchisel and evaluate with compute_cai_validated.

    CAI is always recomputed with compute_cai_validated — DNAchisel's
    own CAI output is NOT trusted.
    """
    if not _DNACHISEL_AVAILABLE:
        return GeneBenchmarkResult(
            gene_name="", organism=organism, tool="dnachisel",
            protein_length=len(protein), cai=0.0, gc_content=0.0,
            runtime_s=0.0, success=False, error="DNAchisel not installed",
        )

    species = _SPECIES_MAP.get(organism, "e_coli")

    # Build initial sequence from preferred codons
    preferred = PREFERRED_CODON_TABLES.get(organism, {})
    if preferred:
        start_seq = "".join(preferred.get(aa, AA_TO_CODONS[aa][0]) for aa in protein)
    else:
        start_seq = "".join(AA_TO_CODONS[aa][0] for aa in protein)

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
        elapsed = time.perf_counter() - t0

        seq = str(problem.sequence)
        # CRITICAL: Recompute CAI with validated evaluator
        cai = compute_cai_validated(seq, organism)
        gc = gc_content(seq)

        # Check constraint violations
        violations = 0
        # GC violation
        if not (GC_LO <= gc <= GC_HI):
            violations += 1

        return GeneBenchmarkResult(
            gene_name="", organism=organism, tool="dnachisel",
            protein_length=len(protein), cai=cai, gc_content=gc,
            runtime_s=elapsed, success=True,
            constraint_violations=violations, sequence=seq,
        )
    except Exception as e:
        elapsed = time.perf_counter() - t0
        return GeneBenchmarkResult(
            gene_name="", organism=organism, tool="dnachisel",
            protein_length=len(protein), cai=0.0, gc_content=0.0,
            runtime_s=elapsed, success=False, error=str(e),
        )


# ═══════════════════════════════════════════════════════════════════════
# Statistical Analysis
# ═══════════════════════════════════════════════════════════════════════

def compute_paired_ttest(bc_values: list[float], dc_values: list[float]) -> tuple[float, float]:
    """Compute paired t-test, return (t_statistic, p_value)."""
    try:
        from scipy import stats as sp_stats
        import numpy as np
        if len(bc_values) < 2:
            return 0.0, 1.0
        bc_arr = np.array(bc_values, dtype=np.float64)
        dc_arr = np.array(dc_values, dtype=np.float64)
        t_stat, p_val = sp_stats.ttest_rel(bc_arr, dc_arr)
        return float(t_stat), float(p_val)
    except Exception:
        return 0.0, 1.0


# ═══════════════════════════════════════════════════════════════════════
# Plotting
# ═══════════════════════════════════════════════════════════════════════

def plot_cai_comparison(
    results: list[GeneBenchmarkResult],
    output_path: Path,
) -> None:
    """Generate CAI comparison bar chart grouped by organism."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import numpy as np

    # Organize paired data by organism
    organisms: dict[str, dict[str, dict[str, float]]] = {}
    for r in results:
        if not r.success:
            continue
        org = r.organism
        if org not in organisms:
            organisms[org] = {}
        organisms[org][r.gene_name] = organisms[org].get(r.gene_name, {})
        organisms[org][r.gene_name][r.tool] = r.cai

    org_list = sorted(organisms.keys())
    n_orgs = len(org_list)

    fig, axes = plt.subplots(1, n_orgs, figsize=(6 * n_orgs, 7), sharey=True)
    fig.suptitle(
        "CAI Comparison: BioCompiler vs DNAchisel\n"
        "(Computed with compute_cai_validated — Sharp & Li 1987)",
        fontsize=14, fontweight="bold",
    )

    colors_bc = "#2ecc71"
    colors_dc = "#e74c3c"

    for idx, org in enumerate(org_list):
        ax = axes[idx] if n_orgs > 1 else axes
        data = organisms[org]

        # Only keep genes with paired data
        genes = sorted(g for g in data if "biocompiler" in data[g] and "dnachisel" in data[g])
        if not genes:
            ax.set_title(f"{org.replace('_', ' ')}\nNo paired data")
            continue

        bc_cais = [data[g]["biocompiler"] for g in genes]
        dc_cais = [data[g]["dnachisel"] for g in genes]

        x = np.arange(len(genes))
        width = 0.35

        ax.bar(x - width/2, bc_cais, width, label="BioCompiler", color=colors_bc, alpha=0.85)
        ax.bar(x + width/2, dc_cais, width, label="DNAchisel", color=colors_dc, alpha=0.85)

        ax.set_xlabel("Gene")
        ax.set_ylabel("CAI")
        ax.set_title(f"{org.replace('_', ' ')}\n(n={len(genes)})")
        ax.set_xticks(x)
        short_names = [g[:10] for g in genes]
        ax.set_xticklabels(short_names, rotation=45, ha="right", fontsize=7)
        ax.legend(fontsize=9)
        ax.set_ylim(0, 1.05)
        ax.grid(axis="y", alpha=0.3)

    plt.tight_layout()
    fig.savefig(output_path, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"  CAI comparison plot saved to {output_path}")


def plot_speed_comparison(
    results: list[GeneBenchmarkResult],
    output_path: Path,
) -> None:
    """Generate speed comparison scatter plot."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import numpy as np

    # Organize paired data
    bc_data: dict[str, tuple[int, float]] = {}
    dc_data: dict[str, tuple[int, float]] = {}
    for r in results:
        if not r.success:
            continue
        key = f"{r.gene_name}_{r.organism}"
        if r.tool == "biocompiler":
            bc_data[key] = (r.protein_length, r.runtime_s * 1000)
        else:
            dc_data[key] = (r.protein_length, r.runtime_s * 1000)

    fig, ax = plt.subplots(figsize=(10, 6))
    fig.suptitle("Runtime Comparison: BioCompiler vs DNAchisel", fontsize=14, fontweight="bold")

    common_keys = sorted(set(bc_data.keys()) & set(dc_data.keys()))

    bc_lengths = [bc_data[k][0] for k in common_keys]
    bc_times = [bc_data[k][1] for k in common_keys]
    dc_lengths = [dc_data[k][0] for k in common_keys]
    dc_times = [dc_data[k][1] for k in common_keys]

    ax.scatter(bc_lengths, bc_times, c="#2ecc71", s=60, alpha=0.7,
               label="BioCompiler", zorder=3, marker="o")
    ax.scatter(dc_lengths, dc_times, c="#e74c3c", s=60, alpha=0.7,
               marker="^", label="DNAchisel", zorder=3)

    # Draw lines connecting paired results
    for k in common_keys:
        bl, bt = bc_data[k]
        dl, dt = dc_data[k]
        ax.plot([bl, dl], [bt, dt], color="gray", alpha=0.3, linewidth=0.5)

    ax.set_xlabel("Protein Length (amino acids)")
    ax.set_ylabel("Runtime (ms)")
    ax.legend()
    ax.grid(alpha=0.3)
    ax.set_yscale("log")

    plt.tight_layout()
    fig.savefig(output_path, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"  Speed comparison plot saved to {output_path}")


def plot_dashboard(
    results: list[GeneBenchmarkResult],
    stats: AggregateStats,
    output_path: Path,
) -> None:
    """Generate a summary dashboard with key metrics."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import numpy as np

    fig, axes = plt.subplots(2, 2, figsize=(14, 10))
    fig.suptitle(
        "BioCompiler vs DNAchisel — Comprehensive E2E Benchmark Dashboard",
        fontsize=15, fontweight="bold",
    )

    # ── Panel 1: Per-gene CAI delta ──
    ax1 = axes[0, 0]
    gene_pairs: dict[str, dict[str, float]] = {}
    for r in results:
        if not r.success:
            continue
        key = f"{r.gene_name}\n({r.organism[:4]})"
        if key not in gene_pairs:
            gene_pairs[key] = {}
        gene_pairs[key][r.tool] = r.cai

    genes_sorted = sorted(gene_pairs.keys())
    deltas = []
    for g in genes_sorted:
        bc = gene_pairs[g].get("biocompiler", 0)
        dc = gene_pairs[g].get("dnachisel", 0)
        deltas.append(bc - dc)

    colors = ["#2ecc71" if d > 0.001 else "#e74c3c" if d < -0.001 else "#95a5a6" for d in deltas]
    ax1.bar(range(len(deltas)), deltas, color=colors, alpha=0.85)
    ax1.axhline(y=0, color="black", linewidth=0.8)
    ax1.set_xticks(range(len(genes_sorted)))
    ax1.set_xticklabels(genes_sorted, rotation=90, fontsize=6)
    ax1.set_ylabel("CAI Delta (BC - DC)")
    ax1.set_title("Per-Gene CAI Advantage (BioCompiler - DNAchisel)")
    ax1.grid(axis="y", alpha=0.3)

    # ── Panel 2: Speed ratio by protein length ──
    ax2 = axes[0, 1]
    bc_dict: dict[str, tuple[int, float]] = {}
    dc_dict: dict[str, tuple[int, float]] = {}
    for r in results:
        if not r.success:
            continue
        key = f"{r.gene_name}_{r.organism}"
        if r.tool == "biocompiler":
            bc_dict[key] = (r.protein_length, r.runtime_s)
        else:
            dc_dict[key] = (r.protein_length, r.runtime_s)

    common = set(bc_dict.keys()) & set(dc_dict.keys())
    lengths = [bc_dict[k][0] for k in common if dc_dict[k][1] > 0]
    ratios = [bc_dict[k][1] / dc_dict[k][1] for k in common if dc_dict[k][1] > 0]

    ax2.scatter(lengths, ratios, c="#3498db", s=60, alpha=0.7)
    ax2.axhline(y=1.0, color="red", linewidth=1, linestyle="--", label="1.0x (equal speed)")
    ax2.set_xlabel("Protein Length (aa)")
    ax2.set_ylabel("Speed Ratio (BC/DC)")
    ax2.set_title("Speed Ratio vs Protein Length")
    ax2.legend()
    ax2.grid(alpha=0.3)

    # ── Panel 3: Win/Loss/Tie pie ──
    ax3 = axes[1, 0]
    win_counts = [stats.bc_wins, stats.dc_wins, stats.ties]
    labels = [f"BC Wins ({stats.bc_wins})", f"DC Wins ({stats.dc_wins})", f"Ties ({stats.ties})"]
    colors_pie = ["#2ecc71", "#e74c3c", "#95a5a6"]
    if sum(win_counts) > 0:
        ax3.pie(win_counts, labels=labels, colors=colors_pie, autopct="%1.1f%%", startangle=90)
    ax3.set_title("Head-to-Head Wins (CAI)")

    # ── Panel 4: Summary stats text ──
    ax4 = axes[1, 1]
    ax4.axis("off")
    summary_text = (
        f"AGGREGATE STATISTICS\n"
        f"{'='*40}\n\n"
        f"Genes tested: {stats.num_genes}\n\n"
        f"Mean CAI:\n"
        f"  BioCompiler: {stats.mean_bc_cai:.4f}\n"
        f"  DNAchisel:   {stats.mean_dc_cai:.4f}\n"
        f"  Delta:       {stats.mean_bc_cai - stats.mean_dc_cai:+.4f}\n\n"
        f"Mean Runtime:\n"
        f"  BioCompiler: {stats.mean_bc_time_ms:.1f} ms\n"
        f"  DNAchisel:   {stats.mean_dc_time_ms:.1f} ms\n"
        f"  Speed ratio: {stats.mean_speed_ratio:.2f}x\n\n"
        f"Head-to-Head:\n"
        f"  BC wins: {stats.bc_wins}  |  DC wins: {stats.dc_wins}  |  Ties: {stats.ties}\n\n"
        f"Paired t-test:\n"
        f"  t = {stats.paired_t_statistic:.3f}\n"
        f"  p = {stats.paired_t_pvalue:.6f}\n"
        f"  {'Significant (p < 0.05)' if stats.paired_t_pvalue < 0.05 else 'Not significant (p >= 0.05)'}"
    )
    ax4.text(0.05, 0.95, summary_text, transform=ax4.transAxes,
             fontsize=11, verticalalignment="top", fontfamily="monospace",
             bbox=dict(boxstyle="round", facecolor="lightyellow", alpha=0.8))

    plt.tight_layout()
    fig.savefig(output_path, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"  Dashboard plot saved to {output_path}")


# ═══════════════════════════════════════════════════════════════════════
# Markdown Report Generation
# ═══════════════════════════════════════════════════════════════════════

def generate_markdown_report(
    results: list[GeneBenchmarkResult],
    stats: AggregateStats,
    output_path: Path,
) -> None:
    """Generate Markdown summary with the required table format."""
    lines = []

    lines.append("# Comprehensive E2E Benchmark: BioCompiler vs DNAchisel")
    lines.append("")
    lines.append("## Fairness Guarantee")
    lines.append("")
    lines.append("Both tools are evaluated with `compute_cai_validated` — the same validated")
    lines.append("CAI function following Sharp & Li (1987). DNAchisel's own CAI output is")
    lines.append("**NOT trusted**. DNAchisel uses `CodonOptimize(species=...)` as an objective,")
    lines.append("but CAI is always recomputed with our validated evaluator.")
    lines.append("")

    # ── Main results table ──
    lines.append("## Per-Gene Results")
    lines.append("")
    lines.append("| Gene | Organism | BC CAI | DC CAI | BC Time (ms) | DC Time (ms) | Speed Ratio | BC Wins |")
    lines.append("|------|----------|--------|--------|-------------|-------------|-------------|---------|")

    # Group by gene for paired display
    gene_pairs: dict[str, dict[str, GeneBenchmarkResult]] = {}
    for r in results:
        key = f"{r.gene_name}_{r.organism}"
        if key not in gene_pairs:
            gene_pairs[key] = {}
        gene_pairs[key][r.tool] = r

    for key in sorted(gene_pairs.keys()):
        pair = gene_pairs[key]
        bc = pair.get("biocompiler")
        dc = pair.get("dnachisel")

        if bc is None or dc is None:
            continue
        if not bc.success and not dc.success:
            continue

        bc_cai = bc.cai if bc.success else 0.0
        dc_cai = dc.cai if dc.success else 0.0
        bc_time = bc.runtime_s * 1000 if bc.success else 0.0
        dc_time = dc.runtime_s * 1000 if dc.success else 0.0

        speed_ratio = bc_time / dc_time if dc_time > 0 else float("inf")
        bc_wins = "[PASS]" if bc_cai > dc_cai + 0.001 else ("tie" if abs(bc_cai - dc_cai) <= 0.001 else "")

        org_display = bc.organism.replace("Escherichia_coli", "E. coli").replace(
            "Homo_sapiens", "Human"
        ).replace("Saccharomyces_cerevisiae", "Yeast")

        lines.append(
            f"| {bc.gene_name} | {org_display} | {bc_cai:.4f} | {dc_cai:.4f} | "
            f"{bc_time:.1f} | {dc_time:.1f} | {speed_ratio:.2f}x | {bc_wins} |"
        )

    lines.append("")

    # ── Aggregate statistics ──
    lines.append("## Aggregate Statistics")
    lines.append("")
    lines.append(f"- **Mean CAI across all genes**: BioCompiler = {stats.mean_bc_cai:.4f}, "
                 f"DNAchisel = {stats.mean_dc_cai:.4f}")
    lines.append(f"- **Mean speed ratio**: {stats.mean_speed_ratio:.2f}x "
                 f"({'BC faster' if stats.mean_speed_ratio < 1.0 else 'DC faster'})")
    lines.append(f"- **Head-to-head wins**: BC = {stats.bc_wins}, DC = {stats.dc_wins}, "
                 f"Ties = {stats.ties}")
    lines.append(f"- **Paired t-test p-value**: {stats.paired_t_pvalue:.6f} "
                 f"({'significant at p<0.05' if stats.paired_t_pvalue < 0.05 else 'not significant'})")
    lines.append("")

    # ── Per-organism breakdown ──
    lines.append("## Per-Organism Breakdown")
    lines.append("")

    for org_key in ["Escherichia_coli", "Homo_sapiens", "Saccharomyces_cerevisiae"]:
        org_results = {k: v for k, v in gene_pairs.items() if org_key in k}
        if not org_results:
            continue

        org_display = org_key.replace("Escherichia_coli", "E. coli").replace(
            "Homo_sapiens", "Human"
        ).replace("Saccharomyces_cerevisiae", "Yeast")

        bc_cais = []
        dc_cais = []
        bc_times = []
        dc_times = []
        bc_wins_count = 0
        dc_wins_count = 0

        for key, pair in org_results.items():
            bc = pair.get("biocompiler")
            dc = pair.get("dnachisel")
            if bc and bc.success:
                bc_cais.append(bc.cai)
                bc_times.append(bc.runtime_s * 1000)
            if dc and dc.success:
                dc_cais.append(dc.cai)
                dc_times.append(dc.runtime_s * 1000)
            if bc and dc and bc.success and dc.success:
                if bc.cai > dc.cai + 0.001:
                    bc_wins_count += 1
                elif dc.cai > bc.cai + 0.001:
                    dc_wins_count += 1

        lines.append(f"### {org_display}")
        lines.append("")
        if bc_cais:
            lines.append(f"- BC mean CAI: {statistics.mean(bc_cais):.4f}")
        if dc_cais:
            lines.append(f"- DC mean CAI: {statistics.mean(dc_cais):.4f}")
        if bc_times:
            lines.append(f"- BC mean time: {statistics.mean(bc_times):.1f} ms")
        if dc_times:
            lines.append(f"- DC mean time: {statistics.mean(dc_times):.1f} ms")
        lines.append(f"- BC wins: {bc_wins_count}, DC wins: {dc_wins_count}")
        lines.append("")

    output_path.write_text("\n".join(lines), encoding="utf-8")
    print(f"  Markdown report saved to {output_path}")


# ═══════════════════════════════════════════════════════════════════════
# Main Benchmark Runner
# ═══════════════════════════════════════════════════════════════════════

def run_comprehensive_benchmark(output_dir: Path) -> None:
    """Run the full comprehensive benchmark."""
    output_dir.mkdir(parents=True, exist_ok=True)

    print("=" * 90)
    print("  COMPREHENSIVE E2E BENCHMARK: BioCompiler vs DNAchisel")
    print("  CAI computed with compute_cai_validated for BOTH tools")
    print("=" * 90)
    print()

    if _DNACHISEL_AVAILABLE:
        print("  DNAchisel: AVAILABLE (CodonOptimize objective + validated CAI)")
    else:
        print("  DNAchisel: NOT INSTALLED — BioCompiler-only results")
    print()

    # ── Run benchmarks ──
    all_results: list[GeneBenchmarkResult] = []

    # Print header
    print(f"  {'Gene':<18} {'Organism':<24} {'Len':>4} │ "
          f"{'BC CAI':>8} {'DC CAI':>8} {'Δ':>8} │ "
          f"{'BC ms':>8} {'DC ms':>8} {'Ratio':>6} │ {'BC Win':>6}")
    print("  " + "-" * 86)

    for entry in PROTEIN_PANEL:
        name = entry["name"]
        protein = entry["protein"]
        organism = entry["organism"]
        aa_len = len(protein)

        # BioCompiler
        bc_result = optimize_biocompiler(protein, organism)
        bc_result.gene_name = name
        all_results.append(bc_result)

        # DNAchisel
        dc_result = optimize_dnachisel(protein, organism)
        dc_result.gene_name = name
        all_results.append(dc_result)

        # Print results
        bc_cai = bc_result.cai if bc_result.success else 0.0
        dc_cai = dc_result.cai if dc_result.success else 0.0
        delta = bc_cai - dc_cai
        bc_ms = bc_result.runtime_s * 1000
        dc_ms = dc_result.runtime_s * 1000
        ratio = bc_ms / dc_ms if dc_ms > 0 else float("inf")
        bc_win = "[PASS]" if delta > 0.001 else ("tie" if abs(delta) <= 0.001 else "")

        org_short = organism.replace("Escherichia_coli", "E.coli").replace(
            "Homo_sapiens", "Human"
        ).replace("Saccharomyces_cerevisiae", "Yeast")

        err_flag = ""
        if not bc_result.success:
            err_flag += " [BC ERR]"
        if not dc_result.success:
            err_flag += " [DC ERR]"

        ratio_str = f"{ratio:.2f}x" if ratio != float("inf") else "N/A"

        print(f"  {name:<18} {org_short:<24} {aa_len:>4} │ "
              f"{bc_cai:>8.4f} {dc_cai:>8.4f} {delta:>+8.4f} │ "
              f"{bc_ms:>8.1f} {dc_ms:>8.1f} {ratio_str:>6} │ {bc_win:>6}{err_flag}")

    print("  " + "-" * 86)
    print()

    # ── Compute aggregate statistics ──
    # Pair results by gene_organism key
    gene_pairs: dict[str, dict[str, GeneBenchmarkResult]] = {}
    for r in all_results:
        key = f"{r.gene_name}_{r.organism}"
        if key not in gene_pairs:
            gene_pairs[key] = {}
        gene_pairs[key][r.tool] = r

    bc_cais_paired = []
    dc_cais_paired = []
    bc_times_paired = []
    dc_times_paired = []
    bc_wins = 0
    dc_wins = 0
    ties = 0

    for key, pair in gene_pairs.items():
        bc = pair.get("biocompiler")
        dc = pair.get("dnachisel")
        if bc and dc and bc.success and dc.success:
            bc_cais_paired.append(bc.cai)
            dc_cais_paired.append(dc.cai)
            bc_times_paired.append(bc.runtime_s)
            dc_times_paired.append(dc.runtime_s)
            if bc.cai > dc.cai + 0.001:
                bc_wins += 1
            elif dc.cai > bc.cai + 0.001:
                dc_wins += 1
            else:
                ties += 1

    n_paired = len(bc_cais_paired)
    mean_bc_cai = statistics.mean(bc_cais_paired) if bc_cais_paired else 0.0
    mean_dc_cai = statistics.mean(dc_cais_paired) if dc_cais_paired else 0.0
    mean_bc_time_ms = statistics.mean([t * 1000 for t in bc_times_paired]) if bc_times_paired else 0.0
    mean_dc_time_ms = statistics.mean([t * 1000 for t in dc_times_paired]) if dc_times_paired else 0.0
    mean_speed_ratio = mean_bc_time_ms / mean_dc_time_ms if mean_dc_time_ms > 0 else 0.0

    # Paired t-test on CAI
    t_stat, t_pval = compute_paired_ttest(bc_cais_paired, dc_cais_paired)

    stats = AggregateStats(
        num_genes=n_paired,
        mean_bc_cai=round(mean_bc_cai, 4),
        mean_dc_cai=round(mean_dc_cai, 4),
        mean_bc_time_ms=round(mean_bc_time_ms, 2),
        mean_dc_time_ms=round(mean_dc_time_ms, 2),
        mean_speed_ratio=round(mean_speed_ratio, 4),
        bc_wins=bc_wins,
        dc_wins=dc_wins,
        ties=ties,
        paired_t_statistic=round(t_stat, 4),
        paired_t_pvalue=round(t_pval, 6),
    )

    # ── Print summary ──
    print()
    print("=" * 90)
    print("  AGGREGATE STATISTICS")
    print("=" * 90)
    print(f"  Paired genes: {n_paired}")
    print()
    print(f"  Mean CAI:   BioCompiler = {mean_bc_cai:.4f}  |  DNAchisel = {mean_dc_cai:.4f}  |  "
          f"Delta = {mean_bc_cai - mean_dc_cai:+.4f}")
    print(f"  Mean Time:  BioCompiler = {mean_bc_time_ms:.1f} ms  |  DNAchisel = {mean_dc_time_ms:.1f} ms  |  "
          f"Ratio = {mean_speed_ratio:.2f}x")
    print()
    print(f"  Head-to-Head:  BC wins = {bc_wins}  |  DC wins = {dc_wins}  |  Ties = {ties}")
    print()
    print(f"  Paired t-test (CAI):  t = {t_stat:.3f}  |  p = {t_pval:.6f}  |  "
          f"{'SIGNIFICANT (p < 0.05)' if t_pval < 0.05 else 'not significant (p >= 0.05)'}")

    # ── Per-organism breakdown ──
    print()
    print("  PER-ORGANISM BREAKDOWN")
    print("  " + "-" * 60)

    for org_key, org_label in [
        ("Escherichia_coli", "E. coli"),
        ("Homo_sapiens", "Human"),
        ("Saccharomyces_cerevisiae", "Yeast"),
    ]:
        org_pairs = {k: v for k, v in gene_pairs.items() if org_key in k}
        if not org_pairs:
            continue

        org_bc_cais = []
        org_dc_cais = []
        org_bc_wins = 0
        org_dc_wins = 0

        for key, pair in org_pairs.items():
            bc = pair.get("biocompiler")
            dc = pair.get("dnachisel")
            if bc and bc.success:
                org_bc_cais.append(bc.cai)
            if dc and dc.success:
                org_dc_cais.append(dc.cai)
            if bc and dc and bc.success and dc.success:
                if bc.cai > dc.cai + 0.001:
                    org_bc_wins += 1
                elif dc.cai > bc.cai + 0.001:
                    org_dc_wins += 1

        bc_mean = statistics.mean(org_bc_cais) if org_bc_cais else 0.0
        dc_mean = statistics.mean(org_dc_cais) if org_dc_cais else 0.0
        print(f"  {org_label:<12}: BC CAI = {bc_mean:.4f}  |  DC CAI = {dc_mean:.4f}  |  "
              f"BC wins = {org_bc_wins}  |  DC wins = {org_dc_wins}")

    # ── Overall verdict ──
    print()
    if bc_wins > dc_wins:
        print(f"  VERDICT: BioCompiler leads with {bc_wins}/{n_paired} head-to-head CAI wins")
    elif dc_wins > bc_wins:
        print(f"  VERDICT: DNAchisel leads with {dc_wins}/{n_paired} head-to-head CAI wins")
    else:
        print(f"  VERDICT: Tie — each tool wins {bc_wins}/{n_paired} comparisons")
    print("=" * 90)
    print()

    # ── Save JSON results ──
    json_data = {
        "metadata": {
            "benchmark_type": "comprehensive_e2e",
            "num_proteins": len(PROTEIN_PANEL),
            "organisms": list(set(e["organism"] for e in PROTEIN_PANEL)),
            "cai_evaluator": "compute_cai_validated (Sharp & Li 1987)",
            "dnachisel_available": _DNACHISEL_AVAILABLE,
            "gc_range": [GC_LO, GC_HI],
        },
        "aggregate_stats": asdict(stats),
        "per_gene_results": [],
    }

    for r in all_results:
        d = {
            "gene_name": r.gene_name,
            "organism": r.organism,
            "tool": r.tool,
            "protein_length": r.protein_length,
            "cai": round(r.cai, 6),
            "gc_content": round(r.gc_content, 6),
            "runtime_s": round(r.runtime_s, 6),
            "success": r.success,
            "constraint_violations": r.constraint_violations,
            "error": r.error,
        }
        json_data["per_gene_results"].append(d)

    json_path = output_dir / "comprehensive_e2e_results.json"
    with open(json_path, "w") as f:
        json.dump(json_data, f, indent=2)
    print(f"  JSON results saved to {json_path}")

    # ── Generate Markdown report ──
    md_path = output_dir / "COMPREHENSIVE_E2E_REPORT.md"
    generate_markdown_report(all_results, stats, md_path)

    # ── Generate plots ──
    try:
        cai_plot_path = output_dir / "comprehensive_cai_comparison.png"
        plot_cai_comparison(all_results, cai_plot_path)
    except Exception as e:
        print(f"  Warning: Could not generate CAI plot: {e}")

    try:
        speed_plot_path = output_dir / "comprehensive_speed_comparison.png"
        plot_speed_comparison(all_results, speed_plot_path)
    except Exception as e:
        print(f"  Warning: Could not generate speed plot: {e}")

    try:
        dashboard_path = output_dir / "comprehensive_dashboard.png"
        plot_dashboard(all_results, stats, dashboard_path)
    except Exception as e:
        print(f"  Warning: Could not generate dashboard plot: {e}")

    print()
    print("  BENCHMARK COMPLETE")
    print()


# ═══════════════════════════════════════════════════════════════════════
# Entry Point
# ═══════════════════════════════════════════════════════════════════════

def main() -> None:
    parser = argparse.ArgumentParser(
        description="Comprehensive E2E benchmark: BioCompiler vs DNAchisel",
    )
    parser.add_argument(
        "--output-dir",
        default="benchmark_results",
        help="Output directory for results (default: benchmark_results)",
    )
    args = parser.parse_args()

    output_dir = Path(args.output_dir)
    if not output_dir.is_absolute():
        output_dir = Path(__file__).parent.parent / output_dir

    run_comprehensive_benchmark(output_dir)


if __name__ == "__main__":
    main()
