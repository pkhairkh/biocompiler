#!/usr/bin/env python3
"""
Comprehensive Head-to-Head Benchmark: BioCompiler vs DNAchisel
================================================================

Fair, statistically rigorous comparison across:
  - 10+ proteins of varying lengths (50aa to 500aa)
  - Multiple organisms (E. coli, Human, Yeast)
  - Both CAI and speed metrics
  - Statistical comparison (paired t-test, Wilcoxon signed-rank)

KEY FAIRNESS GUARANTEE:
  Both tools are evaluated with compute_cai_validated — the same
  validated CAI function.  DNAchisel's own CAI output is NOT trusted.

Outputs:
  - JSON results file (benchmark_results/comprehensive_results.json)
  - Markdown summary with tables (benchmark_results/comprehensive_summary.md)
  - CAI comparison plot (benchmark_results/cai_comparison.png)
  - Speed comparison plot (benchmark_results/speed_comparison.png)

Usage:
    python scripts/comprehensive_head_to_head.py
    python scripts/comprehensive_head_to_head.py --output-dir /tmp/bench
"""

from __future__ import annotations

import argparse
import json
import logging
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
# Protein Panel — 10+ proteins of varying lengths across 3 organisms
# ═══════════════════════════════════════════════════════════════════════

PROTEIN_PANEL: list[dict[str, Any]] = [
    # ── E. coli targets (prokaryotic) ──
    {
        "name": "trpR_Ecoli",
        "protein": "MSIKELIEVQKGIVAKLEQFLPPVEQIKRILRKPGELKRLTLEQYAKQMEEAIRKLTEQFDKDVPQRKRVLLDTTNLLAEIEKELRAQIDNIIEDVDALNKLKKEIEGELQKRLEQVIEKLMDAVAKRRRLSSNAIRKRLADYVK",
        "organism": "Escherichia_coli",
        "length_category": "short",
    },
    {
        "name": "infA_Ecoli",
        "protein": "MASELIRKLAENAIKQAGFPEVMDAFRSQVNELLEKGFQIQVGFPSSKPQTDVDAVKLLEKQGRKRVVAPFIDRGAEKVIKAYSKGVKKPIKPQFADGMLGTVLTGKGRVKAVSGMVGAPGGAKGAKKV",
        "organism": "Escherichia_coli",
        "length_category": "short",
    },
    {
        "name": "recA_Ecoli",
        "protein": "MNLTELKNKPVMEIAQEQLLEQGKIVRALEQGFPDRHVNVLDFIDRAKQIDRGIERDVIVFCDGKERELAAIYQRTLTRGHPQGFGYAQGRFLQGEMAEAWRRADENRDAVFRHLSDPSQPLVVFGMKTTDAVKRIVVDAFQMTRQEAVEAWEKLAQALREGRGIEVIPRSADLIEATRRLFGQHVVNMAAGEGKTVNAVNLGLRGVEVQRLRLDLAFPGELDRAVARLADAREVIHPDVKRLAAGLAEGVRRYLDISVRRMADAIERQEAK",
        "organism": "Escherichia_coli",
        "length_category": "medium",
    },
    {
        "name": "groEL_Ecoli",
        "protein": "MAKDVKFNGELVKFANDAVKVMLEQKPVTVLEQGMKDLRAINILKDAKVKGFKGEVKQIDKLGDGILVSAVGPKTEALVEALKQYVETLADKVGRSVQVLDAVQEFNELEGWKVQGETQLEVKDQIVTKDAFETLDEKGLQKLKNEMQRLDAGKILVTGVGQTEAHVDAKLNRVDMLMDKLVEAGVKVAGTVIDLGKASAEADKLLKELEKGVKETVLPGGVVLTVADKAGLQAEVKEMEKLQDKVKARLEGVVVDTAVPAPVKELVQKMVKEMDQEKLQERIRAALEKAKELVKTRIAEEVKDALKDKAPLVDVKKEIEKRGIESKIIDKVIVAKVAK",
        "organism": "Escherichia_coli",
        "length_category": "medium",
    },
    {
        "name": "tufA_Ecoli",
        "protein": "MAKQFSTKFDALVGVIAAAGKLSEKGKKVILFGVDAKKEDEKIDRLVNEVVKGIYPVTSEDFEYEKEKGKKVFLIPNMFEPVAAKILSEEGRRVGFKVKADVAGVASLDEQRRALRDAVAASIVTIKEGIDRLVSEVTGKIVEGSVAKDIKGKSPEEIERVLKNRIEGVNVIAVGAEGTSDTAKDVLASLVKLGDVVYEVDEAGKTYGEGFLRGVSEAMHSGKAVKIIDKVGGEAIQKEIVADAKLKEKGDVVIPEQGKKVTEAFKRMLQGVDATLFDTDKVIQKVDGNAGDRLAVVCELMDDGKRLPVKVKGYQKLGANERKNKAPQIIVTKYHPDINKNIDVGWLGDKPDAAPLFDLNEKDVRDFVKGKPVTVAVTDGKVNVSAAGAVAVIKDGGVKFN",
        "organism": "Escherichia_coli",
        "length_category": "long",
    },
    # ── Human targets (eukaryotic) ──
    {
        "name": "INS_Human",
        "protein": "MALWMRLLPLLALLALWGPDPAAAFVNQHLCGSHLVEALYLVCGERGFFYTPKTRREAEDLQVGQVELGGGPGAGSLQPLALEGSLQKRGIVEQCCTSICSLYQLENYCN",
        "organism": "Homo_sapiens",
        "length_category": "short",
    },
    {
        "name": "GH1_Human",
        "protein": "MATGSRTSLLLAFGLLCLPWLQEGSAFPTIPLSRLFDNAMLRAHRLHQLAFDTYQEFEEAYIPKEQKYSFLQNPQTSLCFSESIPTPSNREETQQKSNLELLRISLLLIQSWLEPVQFLRSVFANSLVYGASDSNVYDLLKDLEEGIQTLMGRLEDGSPRTGQIFKQTYSKFDTNSHNDDALLKNYGLLYCFRKDMDKVETFLRIVQCRSVEGSCGF",
        "organism": "Homo_sapiens",
        "length_category": "medium",
    },
    {
        "name": "EPO_Human",
        "protein": "MGVHECPAWLWLLLSLLSLPLGLPVLGAPPRLICDSRVLERYLLEAKEAENITTGCAEHCSLNENITVPDTKVNFYAWKRMEVGQQAVEVWQGLALLSEAVLRGQALLVNSSQPWEPLQLHVDKAVSGLRSLTTLLRALGAQKEAISPPDAASAAPLRTITADTFRKLFRVYSNFLRGKLKLYTGEACRTGDR",
        "organism": "Homo_sapiens",
        "length_category": "medium",
    },
    {
        "name": "IFNA2_Human",
        "protein": "MALTFALLVALLVLSCKSSCSVGCDLPQTHSLGSRRTLMLLAQMRRISLFSCLKDRHDFGFPQEEFGNQFQKAETIPVLHEMIQQIFNLFSTKDSSAAWDETLLDKFYTELYQQLNDLEACVIQGVGVTETPLMKEDSILAVRKYFQRITLYLKEKKYSPCAWEVVRAEIMRSFSLSTNLQESLRSKE",
        "organism": "Homo_sapiens",
        "length_category": "medium",
    },
    # ── Yeast targets ──
    {
        "name": "SUC2_Yeast",
        "protein": "MLLQAFLPFLSATVAAYSMLFNSQTQWLTSSNSSSVLLGSTQSVAPFSLFNTTQFSPATATSNLSTIDAYVPQTSINLSLPPSTVSSSNNSVSSNNNPLVNSGNVSLQNLSVTNSNSSLQSTSNGLSLSTVTNSNGQSSLVTSSNQTSGYLSVTSSNTNSQNPSLSSQSILTVSSSNQSISLSVTSSSNQTVS",
        "organism": "Saccharomyces_cerevisiae",
        "length_category": "medium",
    },
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
    # ── Longer proteins ──
    {
        "name": "ampC_Ecoli",
        "protein": "MKTLLLTAVAGAVLAQPSLDSAQREAWLKELKQHPNPNITLKAFSTQNEGKLAEIADALKDKGEVLVVSTQRGMPVVKMNLFSGEKPDMTLFYKNSNPEGLPVFQGTPVDLKNLGKVIVDSFDDVNAVVRMQHGMAFHNFKQETLTNADVEIAHMLKGMVAFKDSQGPTLEILADELNKAKDQIAKLMGSQAKFQDKVVDALHRQLVAGMDAVKDLPAAMLKAGADVIKGTNPLVLDQVAKQLASQADLKILAAKSPLFKLAQKVKMGLKDTVKPGDALSLEKDLLQKAIEKHGDAVIKVVK",
        "organism": "Escherichia_coli",
        "length_category": "long",
    },
    {
        "name": "SERPINA1_Human",
        "protein": "MPSSVSWGILLLAGLCCLVPVSLAEDPQGDAAQKTDTSHHDQDHPTFNKITPNLAEFAFSLYRQLAHQSNSTNIFFSPVSIATAFAMLSLGTKADTHDEILEGLNFNLTEIPEAQIHEGFQELLRTLNQPDSQLQLTTGNGLFLSEGLKLVDKFLEDVKKLYHSEAFTVNFGDTEEAKKQINDYVEKGTQGKIVDLVKELDRDTVFALVNYIFFKGKWERPFEVKDTEEEDFHVDQVTTVKVPMMKRLGMFNIQHCKKLSSWVLLMKYLGNATAIFFLPDEGKLQHLENELTHDIITKFLENEDRRSASLHLPKLSITGTYDLKSVLGQLGITKVFSNGADLSGVTEEAPLKLSKAVHKAVLTIDEKGTEAAGAMFLEAIPMSIPPEVKFNKPFVFLMIEQNTKSPLFMGKVVNPTQK",
        "organism": "Homo_sapiens",
        "length_category": "long",
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
    error: str = ""
    sequence: str = ""


@dataclass
class StatisticalResult:
    """Statistical comparison result for a metric."""
    metric_name: str
    organism: str
    n_pairs: int
    mean_biocompiler: float
    mean_dnachisel: float
    mean_diff: float
    paired_t_statistic: float = 0.0
    paired_t_pvalue: float = 1.0
    wilcoxon_statistic: float = 0.0
    wilcoxon_pvalue: float = 1.0
    cohens_d: float = 0.0
    biocompiler_wins: int = 0
    dnachisel_wins: int = 0
    ties: int = 0


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
        return GeneBenchmarkResult(
            gene_name="", organism=organism, tool="biocompiler",
            protein_length=len(protein), cai=cai, gc_content=gc,
            runtime_s=elapsed, success=True, sequence=seq,
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

    Uses CodonOptimize objective for species-appropriate optimization.
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
        return GeneBenchmarkResult(
            gene_name="", organism=organism, tool="dnachisel",
            protein_length=len(protein), cai=cai, gc_content=gc,
            runtime_s=elapsed, success=True, sequence=seq,
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

def compute_statistics(
    bc_values: list[float],
    dc_values: list[float],
    metric_name: str,
    organism: str,
) -> StatisticalResult:
    """Compute paired t-test, Wilcoxon signed-rank, and effect sizes."""
    from scipy import stats as sp_stats
    import numpy as np

    n = len(bc_values)
    if n < 2:
        return StatisticalResult(
            metric_name=metric_name, organism=organism, n_pairs=n,
            mean_biocompiler=sum(bc_values) / n if n > 0 else 0.0,
            mean_dnachisel=sum(dc_values) / n if n > 0 else 0.0,
            mean_diff=0.0,
        )

    bc_arr = np.array(bc_values, dtype=np.float64)
    dc_arr = np.array(dc_values, dtype=np.float64)
    diffs = bc_arr - dc_arr

    mean_bc = float(np.mean(bc_arr))
    mean_dc = float(np.mean(dc_arr))
    mean_diff = float(np.mean(diffs))

    # Paired t-test
    t_stat, t_pval = sp_stats.ttest_rel(bc_arr, dc_arr)

    # Wilcoxon signed-rank test
    if np.all(diffs == 0):
        w_stat, w_pval = 0.0, 1.0
    else:
        try:
            result = sp_stats.wilcoxon(bc_arr, dc_arr, alternative="two-sided")
            w_stat, w_pval = float(result.statistic), float(result.pvalue)
        except ValueError:
            w_stat, w_pval = 0.0, 1.0

    # Cohen's d
    std_diff = float(np.std(diffs, ddof=1))
    cohens_d = float(np.mean(diffs)) / std_diff if std_diff > 0 else 0.0

    # Win/tie counts
    bc_wins = int(np.sum(bc_arr > dc_arr + 0.001))
    dc_wins = int(np.sum(dc_arr > bc_arr + 0.001))
    ties = n - bc_wins - dc_wins

    return StatisticalResult(
        metric_name=metric_name,
        organism=organism,
        n_pairs=n,
        mean_biocompiler=round(mean_bc, 6),
        mean_dnachisel=round(mean_dc, 6),
        mean_diff=round(mean_diff, 6),
        paired_t_statistic=round(t_stat, 4),
        paired_t_pvalue=round(t_pval, 6),
        wilcoxon_statistic=round(w_stat, 4),
        wilcoxon_pvalue=round(w_pval, 6),
        cohens_d=round(cohens_d, 4),
        biocompiler_wins=bc_wins,
        dnachisel_wins=dc_wins,
        ties=ties,
    )


# ═══════════════════════════════════════════════════════════════════════
# Plotting
# ═══════════════════════════════════════════════════════════════════════

def plot_cai_comparison(
    results: list[GeneBenchmarkResult],
    output_path: Path,
) -> None:
    """Generate CAI comparison plot."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import numpy as np

    # Organize data by organism
    organisms: dict[str, dict[str, list[tuple[str, float]]]] = {}
    for r in results:
        if not r.success:
            continue
        org = r.organism
        if org not in organisms:
            organisms[org] = {"biocompiler": [], "dnachisel": []}
        organisms[org][r.tool].append((r.gene_name, r.cai))

    # Create paired data
    fig, axes = plt.subplots(1, 3, figsize=(18, 6), sharey=True)
    fig.suptitle("CAI Comparison: BioCompiler vs DNAchisel\n(Computed with compute_cai_validated)", fontsize=14, fontweight="bold")

    org_list = sorted(organisms.keys())
    colors_bc = "#2ecc71"  # green
    colors_dc = "#e74c3c"  # red

    for idx, org in enumerate(org_list):
        ax = axes[idx] if len(org_list) > 1 else axes
        data = organisms[org]

        # Find paired genes
        bc_dict = dict(data["biocompiler"])
        dc_dict = dict(data["dnachisel"])
        common_genes = sorted(set(bc_dict.keys()) & set(dc_dict.keys()))

        if not common_genes:
            ax.set_title(f"{org}\nNo paired data")
            continue

        bc_cais = [bc_dict[g] for g in common_genes]
        dc_cais = [dc_dict[g] for g in common_genes]

        x = np.arange(len(common_genes))
        width = 0.35

        bars1 = ax.bar(x - width/2, bc_cais, width, label="BioCompiler", color=colors_bc, alpha=0.8)
        bars2 = ax.bar(x + width/2, dc_cais, width, label="DNAchisel", color=colors_dc, alpha=0.8)

        ax.set_xlabel("Gene")
        ax.set_ylabel("CAI")
        ax.set_title(f"{org}\n(n={len(common_genes)})")
        ax.set_xticks(x)
        short_names = [g[:8] for g in common_genes]
        ax.set_xticklabels(short_names, rotation=45, ha="right", fontsize=8)
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
    """Generate speed comparison plot."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import numpy as np

    # Organize paired data by organism
    organisms: dict[str, dict[str, list[tuple[str, float, int]]]] = {}
    for r in results:
        if not r.success:
            continue
        org = r.organism
        if org not in organisms:
            organisms[org] = {"biocompiler": [], "dnachisel": []}
        organisms[org][r.tool].append((r.gene_name, r.runtime_s, r.protein_length))

    fig, ax = plt.subplots(figsize=(10, 6))
    fig.suptitle("Runtime Comparison: BioCompiler vs DNAchisel", fontsize=14, fontweight="bold")

    bc_lengths = []
    bc_times = []
    dc_lengths = []
    dc_times = []

    for org in organisms:
        bc_dict = {g: (t, l) for g, t, l in organisms[org]["biocompiler"]}
        dc_dict = {g: (t, l) for g, t, l in organisms[org]["dnachisel"]}
        common = set(bc_dict.keys()) & set(dc_dict.keys())
        for g in common:
            bt, bl = bc_dict[g]
            dt, dl = dc_dict[g]
            bc_lengths.append(bl)
            bc_times.append(bt * 1000)  # ms
            dc_lengths.append(dl)
            dc_times.append(dt * 1000)  # ms

    ax.scatter(bc_lengths, bc_times, c="#2ecc71", s=60, alpha=0.7, label="BioCompiler", zorder=3)
    ax.scatter(dc_lengths, dc_times, c="#e74c3c", s=60, alpha=0.7, marker="^", label="DNAchisel", zorder=3)

    ax.set_xlabel("Protein Length (amino acids)")
    ax.set_ylabel("Runtime (ms)")
    ax.legend()
    ax.grid(alpha=0.3)
    ax.set_yscale("log")

    plt.tight_layout()
    fig.savefig(output_path, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"  Speed comparison plot saved to {output_path}")


# ═══════════════════════════════════════════════════════════════════════
# Markdown Report Generation
# ═══════════════════════════════════════════════════════════════════════

def generate_markdown_report(
    results: list[GeneBenchmarkResult],
    stats: list[StatisticalResult],
    output_path: Path,
) -> None:
    """Generate Markdown summary with tables."""
    lines = []

    lines.append("# Comprehensive Head-to-Head Benchmark: BioCompiler vs DNAchisel")
    lines.append("")
    lines.append("## Fairness Guarantee")
    lines.append("")
    lines.append("Both tools are evaluated with `compute_cai_validated` — the same validated")
    lines.append("CAI function following Sharp & Li (1987). DNAchisel's own CAI output is")
    lines.append("**NOT trusted**. DNAchisel uses `CodonOptimize(species=...)` as an objective,")
    lines.append("but CAI is always recomputed with our validated evaluator.")
    lines.append("")

    # Per-gene results table
    lines.append("## Per-Gene Results")
    lines.append("")
    lines.append("| Gene | Organism | Length | BC CAI | DC CAI | CAI Δ | BC Time (ms) | DC Time (ms) | BC GC | DC GC |")
    lines.append("|------|----------|--------|--------|--------|-------|-------------|-------------|-------|-------|")

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
        if not bc.success or not dc.success:
            continue

        cai_delta = bc.cai - dc.cai
        lines.append(
            f"| {bc.gene_name} | {bc.organism.replace('_', ' ')} | {bc.protein_length}aa | "
            f"{bc.cai:.4f} | {dc.cai:.4f} | {cai_delta:+.4f} | "
            f"{bc.runtime_s*1000:.1f} | {dc.runtime_s*1000:.1f} | "
            f"{bc.gc_content:.3f} | {dc.gc_content:.3f} |"
        )

    lines.append("")

    # Statistical comparison tables per organism
    lines.append("## Statistical Comparison")
    lines.append("")

    for org in sorted(set(s.organism for s in stats)):
        org_stats = [s for s in stats if s.organism == org]
        lines.append(f"### {org.replace('_', ' ')}")
        lines.append("")
        lines.append("| Metric | n | BC Mean | DC Mean | Mean Δ | t-stat | t p-value | W-stat | W p-value | Cohen's d | BC Wins | DC Wins | Ties |")
        lines.append("|--------|---|---------|---------|--------|--------|-----------|--------|-----------|-----------|---------|---------|------|")

        for s in org_stats:
            lines.append(
                f"| {s.metric_name} | {s.n_pairs} | {s.mean_biocompiler:.4f} | {s.mean_dnachisel:.4f} | "
                f"{s.mean_diff:+.4f} | {s.paired_t_statistic:.3f} | {s.paired_t_pvalue:.4f} | "
                f"{s.wilcoxon_statistic:.3f} | {s.wilcoxon_pvalue:.4f} | {s.cohens_d:.3f} | "
                f"{s.biocompiler_wins} | {s.dnachisel_wins} | {s.ties} |"
            )
        lines.append("")

    # Overall summary
    lines.append("## Overall Summary")
    lines.append("")
    cai_stats = [s for s in stats if s.metric_name == "CAI"]
    speed_stats = [s for s in stats if s.metric_name == "Runtime"]

    for s_list, label in [(cai_stats, "CAI"), (speed_stats, "Runtime")]:
        total_bc_wins = sum(s.biocompiler_wins for s in s_list)
        total_dc_wins = sum(s.dnachisel_wins for s in s_list)
        total_ties = sum(s.ties for s in s_list)
        total_n = sum(s.n_pairs for s in s_list)
        lines.append(f"**{label}**: BioCompiler wins {total_bc_wins}/{total_n}, "
                     f"DNAchisel wins {total_dc_wins}/{total_n}, Ties {total_ties}/{total_n}")
        lines.append("")

    # Significance summary
    lines.append("### Statistical Significance")
    lines.append("")
    for s in stats:
        sig_marker = ""
        if s.paired_t_pvalue < 0.001:
            sig_marker = " (***)"
        elif s.paired_t_pvalue < 0.01:
            sig_marker = " (**)"
        elif s.paired_t_pvalue < 0.05:
            sig_marker = " (*)"
        effect = "negligible" if abs(s.cohens_d) < 0.2 else "small" if abs(s.cohens_d) < 0.5 else "medium" if abs(s.cohens_d) < 0.8 else "large"
        lines.append(f"- {s.organism.replace('_', ' ')} / {s.metric_name}: "
                     f"Δ={s.mean_diff:+.4f}, t p={s.paired_t_pvalue:.4f}{sig_marker}, "
                     f"W p={s.wilcoxon_pvalue:.4f}, Cohen's d={s.cohens_d:.2f} ({effect})")
    lines.append("")

    output_path.write_text("\n".join(lines), encoding="utf-8")
    print(f"  Markdown report saved to {output_path}")


# ═══════════════════════════════════════════════════════════════════════
# Main Benchmark Runner
# ═══════════════════════════════════════════════════════════════════════

def run_comprehensive_benchmark(output_dir: Path) -> None:
    """Run the full comprehensive benchmark."""
    output_dir.mkdir(parents=True, exist_ok=True)

    print("=" * 80)
    print("COMPREHENSIVE HEAD-TO-HEAD BENCHMARK: BioCompiler vs DNAchisel")
    print("CAI computed with compute_cai_validated for BOTH tools")
    print("=" * 80)
    print()

    if _DNACHISEL_AVAILABLE:
        print("DNAchisel: AVAILABLE (using CodonOptimize objective + validated CAI)")
    else:
        print("DNAchisel: NOT INSTALLED — BioCompiler-only results")
    print()

    # ── Run benchmarks ──
    all_results: list[GeneBenchmarkResult] = []

    print(f"{'Gene':<20} {'Organism':<25} {'Len':>4} │ "
          f"{'BC CAI':>8} {'DC CAI':>8} {'Δ':>8} │ "
          f"{'BC ms':>8} {'DC ms':>8}")
    print("-" * 100)

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

        err_flag = ""
        if not bc_result.success:
            err_flag += " [BC ERR]"
        if not dc_result.success:
            err_flag += " [DC ERR]"

        print(f"{name:<20} {organism.replace('_', ' '):<25} {aa_len:>4} │ "
              f"{bc_cai:>8.4f} {dc_cai:>8.4f} {delta:>+8.4f} │ "
              f"{bc_ms:>8.1f} {dc_ms:>8.1f}{err_flag}")

    print("-" * 100)
    print()

    # ── Statistical analysis per organism ──
    all_stats: list[StatisticalResult] = []

    organisms_with_data = sorted(set(
        r.organism for r in all_results if r.success
    ))

    for org in organisms_with_data:
        # Pair results
        bc_results = {r.gene_name: r for r in all_results
                      if r.organism == org and r.tool == "biocompiler" and r.success}
        dc_results = {r.gene_name: r for r in all_results
                      if r.organism == org and r.tool == "dnachisel" and r.success}

        common_genes = sorted(set(bc_results.keys()) & set(dc_results.keys()))

        if len(common_genes) < 2:
            print(f"  {org}: Not enough paired data for statistics (n={len(common_genes)})")
            continue

        # CAI statistics
        bc_cais = [bc_results[g].cai for g in common_genes]
        dc_cais = [dc_results[g].cai for g in common_genes]
        cai_stats = compute_statistics(bc_cais, dc_cais, "CAI", org)
        all_stats.append(cai_stats)

        # Runtime statistics
        bc_times = [bc_results[g].runtime_s for g in common_genes]
        dc_times = [dc_results[g].runtime_s for g in common_genes]
        time_stats = compute_statistics(bc_times, dc_times, "Runtime", org)
        all_stats.append(time_stats)

        # Print per-organism summary
        print(f"\n  {org.replace('_', ' ')} (n={len(common_genes)}):")
        print(f"    CAI:    BC={cai_stats.mean_biocompiler:.4f}  DC={cai_stats.mean_dnachisel:.4f}  "
              f"Δ={cai_stats.mean_diff:+.4f}  t p={cai_stats.paired_t_pvalue:.4f}  "
              f"W p={cai_stats.wilcoxon_pvalue:.4f}  d={cai_stats.cohens_d:.2f}")
        print(f"    Speed:  BC={time_stats.mean_biocompiler*1000:.1f}ms  DC={time_stats.mean_dnachisel*1000:.1f}ms  "
              f"Δ={time_stats.mean_diff*1000:+.1f}ms  t p={time_stats.paired_t_pvalue:.4f}")

    # ── Save JSON results ──
    json_data = []
    for r in all_results:
        d = asdict(r)
        d["cai"] = round(d["cai"], 6)
        d["gc_content"] = round(d["gc_content"], 6)
        d["runtime_s"] = round(d["runtime_s"], 6)
        json_data.append(d)

    json_path = output_dir / "comprehensive_results.json"
    with open(json_path, "w") as f:
        json.dump({"results": json_data, "statistics": [asdict(s) for s in all_stats]}, f, indent=2)
    print(f"\n  JSON results saved to {json_path}")

    # ── Generate Markdown report ──
    md_path = output_dir / "comprehensive_summary.md"
    generate_markdown_report(all_results, all_stats, md_path)

    # ── Generate plots ──
    try:
        cai_plot_path = output_dir / "cai_comparison.png"
        plot_cai_comparison(all_results, cai_plot_path)
    except Exception as e:
        print(f"  Warning: Could not generate CAI plot: {e}")

    try:
        speed_plot_path = output_dir / "speed_comparison.png"
        plot_speed_comparison(all_results, speed_plot_path)
    except Exception as e:
        print(f"  Warning: Could not generate speed plot: {e}")

    print("\n" + "=" * 80)
    print("BENCHMARK COMPLETE")
    print("=" * 80)


# ═══════════════════════════════════════════════════════════════════════
# Entry Point
# ═══════════════════════════════════════════════════════════════════════

def main() -> None:
    parser = argparse.ArgumentParser(
        description="Comprehensive head-to-head benchmark: BioCompiler vs DNAchisel",
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
