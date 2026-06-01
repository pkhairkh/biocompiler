"""
BioCompiler Comprehensive Benchmark Suite v2.0
================================================

12-gene panel with statistical analysis, ablation study, and Pareto frontier.

Provides:
  1. Multi-gene comparison (12 proteins, both cai_first & constraint_first)
  2. Baselines (SimpleCAI, Random)
  3. Statistical analysis (mean/std/min/max + formatted table)
  4. Ablation study (phase contribution)
  5. Pareto frontier plot (CAI vs constraint violations)
  6. Ablation bar chart

Usage:
  cd /home/z/my-project/biocompiler && \\
  PYTHONPATH="src:$PYTHONPATH" python3 -m biocompiler.comprehensive_benchmark

Output saved to: /home/z/my-project/download/benchmark_results/
"""

import csv
import json
import logging
import math
import os
import random
import time
from datetime import datetime, timezone
from pathlib import Path
from typing import Optional

import matplotlib
matplotlib.use("Agg")
import matplotlib.font_manager as fm
try:
    fm.fontManager.addfont('/usr/share/fonts/truetype/chinese/NotoSansSC[wght].ttf')
except RuntimeError:
    try:
        fm.fontManager.addfont('/usr/share/fonts/truetype/chinese/SarasaMonoSC-Regular.ttf')
    except Exception:
        pass
import matplotlib.pyplot as plt
plt.rcParams['font.sans-serif'] = ['Noto Sans SC', 'DejaVu Sans']
plt.rcParams['axes.unicode_minus'] = False
import numpy as np

from .optimization import BioOptimizer
from .translation import translate, compute_cai
from .scanner import gc_content
from .species import SPECIES
from .constants import CODON_TABLE, AA_TO_CODONS, RESTRICTION_ENZYMES, reverse_complement
from .restriction_sites import get_recognition_site
from .certificates import compute_certificate

logger = logging.getLogger(__name__)

# ============================================================================
# Output Directory
# ============================================================================

OUTPUT_DIR = Path("/home/z/my-project/download/benchmark_results")

# ============================================================================
# 12-Gene Panel (amino acid sequence, organism)
# ============================================================================

GENE_PANEL = {
    "HBB": ("MVHLTPEEKSAVTALWGKVNVDEVGGEALGRLLVVYPWTQRFFESFGDLSTPDAVMGNPKVKAHGKKVLGAFSDGLAHLDNLKGTFATLSELHCDKLHVDPENFRLLGNVLVCVLAHHFGKEFTPPVQAAYQKVVAGVANALAHKYH", "Homo_sapiens"),
    "INS": ("MALWMRLLPLLALLALWGPDPAAAFVNQHLCGSHLVEALYLVCGERGFFYTPKTRREAEDLQVGQVELGGGPGAGSLQPLALEGSLQKRGIVEQCCTSICSLYQLENYCN", "Homo_sapiens"),
    "EGFP": ("MSKGEELFTGVVPILVELDGDVNGHKFSVSGEGEGDATYGKLTLKFICTTGKLPVPWPTLVTTFSYGVQCFSRYPDHMKQHDFFKSAMPEGYVQERTIFFKDDGNYKTRAEVKFEGDTLVNRIELKGIDFKEDGNILGHKLEYNYNSHNVYIMADKQKNGIKVNFKIRHNIEDGSVQLADHYQQNTPIGDGPVLLPDNHYLSTQSALSKDPNEKRDHMVLLEFVTAAGITHGMDELYK", "Homo_sapiens"),
    "TNF_alpha": ("MSTESMIRDVELAEEALPKKTGGPQGSRRCLFLSLFSFLIVAGATTLFCLLNFGVIGPQRKRRPYEIHEVQGVFNITLSCWNYKSSSFSQYLFSRLHDDQNQQIFLKNCSKNSVTWCENLTKSCNIKFNSQICNGRGFCRFHVCSSKGYSRGTIYESESNISKTSYLFQMIQKTSFNSYIFWLHNIKTYNKT", "Homo_sapiens"),
    "IL2": ("MYRMQLLLLSCIALSLALVTNSAPTSSSTKKTQQLELESPSPSPSQDETQLLEHNQLPLSELQELQALQNAVSQSRNLQLESQATLKSLQELQELSQLQKASQVLGQESSFSSYPKLAFSESSKKPSSSQSSSSQFSQSSQFSVQDVVPKLQYQNDVFYFRSKQQYVSNHYSQKTSISP", "Homo_sapiens"),
    "EPO": ("MGVHECPAWLWLLLSLLSLPLGLPVLGAPPRLICDSRVLERYLLEAKEAENITTGCAEHCSLNENITVPDTKVNFYAWKRMEVGQQAVEVWQGLALLSEAVLRGQALLVNSSQPWEPLQLHVDKAVSGLRSLTTLLRALGAQKEAISPPDAASAAPLRTITADTFRKLFRVYSNFLRGKLKLYTGEACRTGDR", "Homo_sapiens"),
    "mCherry": ("MVSKGEEDNMAIIHMFMRFHVMEIESGGDFTYMKKVLYKDNGHIITVEYPNDGKLVEFKFPGDGTIEREHDLFKLEKNKTYLQMLDGMILYVTSGTCLKEDNVKLYKCFHEGIKDANRDLFNDVVTKDTYKLILKVDKHDPSYWKTYQEHPSLFCVKSHPQ", "Homo_sapiens"),
    "IFN_alpha2": ("CDLPQTHSLGNRRTLMLLAQMRKISLFSCLKDRHDFGFPQEEFGNQFQKAETIPVLHEMIQQIFNLFSTKDSSAAWDETLLDKFYTELYQQLNDLEACVIQEVGVQETPLMNEDSILAVRKYFQRITLYLKEKKYSPCAWEVVRAEIMRSFSLSTNLQESLRSKE", "Homo_sapiens"),
    "GH1": ("MPTIPLSRLFDNAMLRAGIVHFCIDKLTNNSSSFSRLFLQGFLNFYSFLQPNGAVFMDSGRQQLLQDYKKKETFYLMKDLEDPQLLRSVLSQDMQHVFYSLLSFQDVFHFVDSCDLVQNYRLSLVSTSMARLRHLVQEYFNLITSFCRKVDHHHMHQNLPQLFQTSRPQPIFSRPILFQKSFTSMLFQNSYQQPQASFPQQPQSQSFPQQPQSQSFPQQPQSQSFPQQPQSQSF", "Homo_sapiens"),
    "BSA_frag": ("MKWVTFISLLLLFSSAYSRGVFRRDTHKSEIAHRFKDLGEEHFKGLVLIAFSQYLQQCPFDEHVKLVNELTEFAKTCVADESHAGCEKSLHTLFGDELCKVASLRETYGDMADCCEKQEPERNECFLSHKDDSPDLPKLKPDPNTLCDEFKADEKKFWGKYLYEIARRHPYFYAPELLYYANKYNGVFQECCQAEDKGACLLPKIETMREKVLTSARQRLRCASIQKFGERALKAWSVARLSQKFPKAEFVEVTKLVTDLTKVHKECCHGDLLECADDRADLAKYICDNQDTISSKLKECCDKPLLEKSHCIAEVEKDAIPENLPPLTADFAEDKDVCKNYQEAKDAFLGSFLYEYSRRHPEYAVSVLLRLAKEYEATLEECCAKDDPHACYSTVFDKLKHLVDEPQNLIKQNCDQFEKLGEYGFQNALIVRYTRKVPQVSTPTLVEVSRSLGKVGTRCCTKPESERMPCTEDYLSLILNRLCVLHEKTPVSEKVTKCCTESLVNRRPCFSALTPDETYVPKAFDEKLFTFHADICTLPDTEKQIKKQTALVELLKHKPKATEEQLKTVMENFVAFVDKCCAADDKEACFAVEGPKLVVSTQTALA", "Homo_sapiens"),
    "LacZ_frag": ("MTMITDSLAVVLQRRDWENPGVTQLNRLAAHPPFASWRNSEEARTDRPSQQLRSLNGEWRFAWFPAPEAVPESWLECDLPEADTVVVPSNWQMHGYDAPIYTNVTYPITVNPPFVPTENPTGCYSLTFNVDESWLQEGQTRIIFDGVNSAFHLWCNGRWVGYGQDSRLPSEFDLSAFLRAGENRLAVMVLRWSDGSYLEDQDMWRMSGIFRDVSLLHKPTTQISDFHVATRFNDDFSRAVLEAEVQMCGELRDYLRVTVSLWQGETQVASGTAPFGGEIIDERGGYADRVTLRLNVENPKLWSAEIPNLYRAVVELHTADGTLIEAEACDVGFREVRIENGLLLLNGKPLLIRGVNRHEHHLGCGSTFDNGSFWTQVRGELGMVDAYRQTRSEGCQIRVQVKVASLPEEATLVLTNDSVFHADAQGWFHPWLSQYF", "Escherichia_coli"),
    "Cas9_frag": ("MDKKYSIGLDIGTNSVGWAVITDEYKVPSKKFKVLGNTDRHSIKKNLIGALLFDSGETAEATRLKRTARRRYTRRKNRICYLQEIFSNEMAKVDDSFFHRLEESFLVEEDKKHERHPIFGNIVDEVAYHEKYPTIYHLRKKLVDSTDKADLRLIYLALAHMIKFRGHFLIEGDLNPDNSDVDKLFIQDVQTGGILKDSKIPAIIRPIFKRKLLFDVYRKNHKAEREKVRMSLDGLIEKFSVKETLKELKKSVIKDNKTIKEVGRRAVNIKKITHVPVEEIARKFDNPMVIKTLEEVKKEEKPVQKIIKKIEEVK", "Streptococcus_pyogenes"),
}

# Mapping from full organism name to BioOptimizer species key
ORGANISM_TO_SPECIES = {
    "Homo_sapiens": "human",
    "Escherichia_coli": "ecoli",
    "Streptococcus_pyogenes": "human",  # fallback: no specific table
}

# Mapping from full organism name to CAI computation organism
# (compute_cai uses the organisms module which has limited support)
ORGANISM_FOR_CAI = {
    "Homo_sapiens": "Homo_sapiens",
    "Escherichia_coli": "Escherichia_coli",
    "Streptococcus_pyogenes": "Homo_sapiens",  # fallback: use human CAI table
}

# Enzymes to check
DEFAULT_ENZYMES = ["EcoRI", "BamHI", "HindIII", "XhoI"]

# Total constraints tracked
MAX_CONSTRAINTS = 5


# ============================================================================
# Helper Functions
# ============================================================================

def _build_best_codon_sequence(protein: str, species: str = "human") -> str:
    """Build initial DNA sequence using highest-CAI codons."""
    usage = SPECIES.get(species, SPECIES["ecoli"])
    result = []
    for aa in protein:
        codons = AA_TO_CODONS.get(aa, [])
        if not codons:
            result.append("NNN")
            continue
        result.append(max(codons, key=lambda c: usage.get(c, 0.0)))
    return "".join(result)


def _count_restriction_sites(sequence: str, enzymes: list[str] | None = None) -> int:
    """Count restriction enzyme recognition sites (both strands)."""
    if not enzymes:
        enzymes = DEFAULT_ENZYMES
    count = 0
    seq_upper = sequence.upper()
    for enz_name in enzymes:
        site = get_recognition_site(enz_name)
        if site is None:
            continue
        site_upper = site.upper()
        if any(b not in "ACGT" for b in site_upper):
            continue
        start = 0
        while True:
            pos = seq_upper.find(site_upper, start)
            if pos == -1:
                break
            count += 1
            start = pos + 1
        site_rc = reverse_complement(site_upper)
        if site_rc != site_upper:
            start = 0
            while True:
                pos = seq_upper.find(site_rc, start)
                if pos == -1:
                    break
                count += 1
                start = pos + 1
    return count


def _count_gt_dinucleotides(seq: str) -> int:
    """Count GT dinucleotides."""
    return sum(1 for i in range(len(seq) - 1) if seq[i:i+2] == "GT")


def _count_cpg_obs_exp(seq: str) -> float:
    """Compute CpG Obs/Exp ratio."""
    c = seq.count("C")
    g = seq.count("G")
    cg = sum(1 for i in range(len(seq) - 1) if seq[i:i+2] == "CG")
    expected = (c * g) / len(seq) if len(seq) > 0 else 0
    return cg / expected if expected > 0 else 0.0


def _satisfied_constraints(
    sequence: str,
    protein: str,
    enzymes: list[str] | None = None,
    gc_lo: float = 0.30,
    gc_hi: float = 0.70,
) -> int:
    """Count satisfied constraints (translation, GC, RS, GT, CpG)."""
    count = 0
    translated = translate(sequence, to_stop=True)
    if translated == protein:
        count += 1
    gc = gc_content(sequence)
    if gc_lo <= gc <= gc_hi:
        count += 1
    if _count_restriction_sites(sequence, enzymes) == 0:
        count += 1
    if _count_gt_dinucleotides(sequence) == 0:
        count += 1
    if _count_cpg_obs_exp(sequence) < 0.6:
        count += 1
    return count


def _compute_metrics(
    sequence: str,
    protein: str,
    organism: str,
    species: str,
    enzymes: list[str],
    gc_lo: float = 0.30,
    gc_hi: float = 0.70,
) -> dict:
    """Compute all metrics for an optimized sequence."""
    cai_organism = ORGANISM_FOR_CAI.get(organism, "Homo_sapiens")
    try:
        cai = compute_cai(sequence, cai_organism)
    except Exception:
        cai = 0.0
    gc = gc_content(sequence)
    rs_count = _count_restriction_sites(sequence, enzymes)
    gt_count = _count_gt_dinucleotides(sequence)
    cpg_ratio = _count_cpg_obs_exp(sequence)
    constraints = _satisfied_constraints(sequence, protein, enzymes, gc_lo, gc_hi)
    violations = MAX_CONSTRAINTS - constraints
    return {
        "cai": cai,
        "gc_content": gc,
        "restriction_site_count": rs_count,
        "gt_count": gt_count,
        "cpg_ratio": round(cpg_ratio, 4),
        "constraints_satisfied": constraints,
        "constraint_violations": violations,
        "max_constraints": MAX_CONSTRAINTS,
    }


# ============================================================================
# Tool Implementations
# ============================================================================

def optimize_biocompiler(
    protein: str,
    organism: str,
    strategy: str = "constraint_first",
    enzymes: list[str] | None = None,
    gc_lo: float = 0.30,
    gc_hi: float = 0.70,
) -> dict:
    """Run BioCompiler with specified strategy (cai_first or constraint_first)."""
    enzymes = enzymes or DEFAULT_ENZYMES
    species = ORGANISM_TO_SPECIES.get(organism, "human")
    tool_name = f"BC_{strategy}"
    t0 = time.perf_counter()
    try:
        opt = BioOptimizer(
            species=species,
            enzymes=enzymes,
            avoid_gt=True,
            strategy=strategy,
        )
        initial_seq = _build_best_codon_sequence(protein, species)
        optimized, pred_results, cert_text = opt.optimize(initial_seq)
        elapsed = time.perf_counter() - t0

        metrics = _compute_metrics(optimized, protein, organism, species, enzymes, gc_lo, gc_hi)
        cert_level = compute_certificate(pred_results)

        return {
            "tool": tool_name,
            "sequence": optimized,
            "execution_time_s": round(elapsed, 4),
            "success": True,
            "certificate_level": cert_level.value if hasattr(cert_level, 'value') else str(cert_level),
            **metrics,
        }
    except Exception as exc:
        elapsed = time.perf_counter() - t0
        logger.error("%s failed: %s", tool_name, exc)
        return {
            "tool": tool_name,
            "sequence": "",
            "execution_time_s": round(elapsed, 4),
            "success": False,
            "error": str(exc),
            "cai": 0.0, "gc_content": 0.0,
            "restriction_site_count": -1, "gt_count": -1,
            "cpg_ratio": 0.0,
            "constraints_satisfied": 0, "constraint_violations": MAX_CONSTRAINTS,
            "max_constraints": MAX_CONSTRAINTS,
        }


def optimize_simple_cai(
    protein: str,
    organism: str,
    enzymes: list[str] | None = None,
) -> dict:
    """SimpleCAI baseline: most-preferred codon only. No constraint handling."""
    enzymes = enzymes or DEFAULT_ENZYMES
    species = ORGANISM_TO_SPECIES.get(organism, "human")
    t0 = time.perf_counter()
    try:
        sequence = _build_best_codon_sequence(protein, species)
        elapsed = time.perf_counter() - t0
        metrics = _compute_metrics(sequence, protein, organism, species, enzymes)
        return {
            "tool": "SimpleCAI",
            "sequence": sequence,
            "execution_time_s": round(elapsed, 4),
            "success": True,
            **metrics,
        }
    except Exception as exc:
        elapsed = time.perf_counter() - t0
        return {
            "tool": "SimpleCAI",
            "sequence": "",
            "execution_time_s": round(elapsed, 4),
            "success": False,
            "error": str(exc),
            "cai": 0.0, "gc_content": 0.0,
            "restriction_site_count": -1, "gt_count": -1,
            "cpg_ratio": 0.0,
            "constraints_satisfied": 0, "constraint_violations": MAX_CONSTRAINTS,
            "max_constraints": MAX_CONSTRAINTS,
        }


def optimize_random(
    protein: str,
    organism: str,
    enzymes: list[str] | None = None,
    seed: int = 42,
) -> dict:
    """Random baseline: frequency-weighted random codon selection."""
    enzymes = enzymes or DEFAULT_ENZYMES
    species = ORGANISM_TO_SPECIES.get(organism, "human")
    t0 = time.perf_counter()
    try:
        rng = random.Random(seed)
        usage = SPECIES.get(species, SPECIES["ecoli"])
        seq_chars = []
        for aa in protein:
            codons = AA_TO_CODONS.get(aa, [])
            if not codons:
                seq_chars.append("NNN")
                continue
            weights = [usage.get(c, 0.01) for c in codons]
            total = sum(weights)
            if total <= 0:
                chosen = rng.choice(codons)
            else:
                probs = [w / total for w in weights]
                chosen = rng.choices(codons, weights=probs, k=1)[0]
            seq_chars.append(chosen)
        sequence = "".join(seq_chars)
        elapsed = time.perf_counter() - t0
        metrics = _compute_metrics(sequence, protein, organism, species, enzymes)
        return {
            "tool": "Random",
            "sequence": sequence,
            "execution_time_s": round(elapsed, 4),
            "success": True,
            **metrics,
        }
    except Exception as exc:
        elapsed = time.perf_counter() - t0
        return {
            "tool": "Random",
            "sequence": "",
            "execution_time_s": round(elapsed, 4),
            "success": False,
            "error": str(exc),
            "cai": 0.0, "gc_content": 0.0,
            "restriction_site_count": -1, "gt_count": -1,
            "cpg_ratio": 0.0,
            "constraints_satisfied": 0, "constraint_violations": MAX_CONSTRAINTS,
            "max_constraints": MAX_CONSTRAINTS,
        }


# ============================================================================
# Part 1: Multi-Gene Comparison
# ============================================================================

def run_multi_gene_comparison(
    genes: dict | None = None,
    enzymes: list[str] | None = None,
    gc_lo: float = 0.30,
    gc_hi: float = 0.70,
) -> list[dict]:
    """Run all tools on all 12 genes."""
    genes = genes or GENE_PANEL
    enzymes = enzymes or DEFAULT_ENZYMES
    all_results = []

    for gene_name, (protein, organism) in genes.items():
        print(f"  Benchmarking {gene_name} ({len(protein)} aa, {organism})...")

        # BioCompiler constraint_first
        bc_cf = optimize_biocompiler(protein, organism, "constraint_first", enzymes, gc_lo, gc_hi)
        bc_cf["gene"] = gene_name
        bc_cf["protein_length"] = len(protein)
        bc_cf["organism"] = organism
        all_results.append(bc_cf)

        # BioCompiler cai_first
        bc_cai = optimize_biocompiler(protein, organism, "cai_first", enzymes, gc_lo, gc_hi)
        bc_cai["gene"] = gene_name
        bc_cai["protein_length"] = len(protein)
        bc_cai["organism"] = organism
        all_results.append(bc_cai)

        # SimpleCAI
        sc = optimize_simple_cai(protein, organism, enzymes)
        sc["gene"] = gene_name
        sc["protein_length"] = len(protein)
        sc["organism"] = organism
        all_results.append(sc)

        # Random
        rn = optimize_random(protein, organism, enzymes)
        rn["gene"] = gene_name
        rn["protein_length"] = len(protein)
        rn["organism"] = organism
        all_results.append(rn)

    return all_results


# ============================================================================
# Part 2: Statistical Analysis
# ============================================================================

def compute_statistics(results: list[dict]) -> dict:
    """Compute descriptive stats: mean, std, min, max for each metric per tool."""
    # Organize by tool
    tools: dict[str, dict[str, list]] = {}
    for r in results:
        tool = r["tool"]
        if tool not in tools:
            tools[tool] = {
                "cai": [], "gc_content": [], "restriction_site_count": [],
                "gt_count": [], "cpg_ratio": [],
                "constraints_satisfied": [], "constraint_violations": [],
                "execution_time_s": [],
            }
        if r.get("success"):
            for metric in tools[tool]:
                val = r.get(metric)
                if val is not None and val >= 0:
                    tools[tool][metric].append(val)

    stats_report = {}
    for tool_name, data in tools.items():
        tool_stats = {}
        for metric_name, values in data.items():
            if not values:
                tool_stats[metric_name] = {"n": 0}
                continue
            arr = np.array(values, dtype=float)
            tool_stats[metric_name] = {
                "n": int(len(arr)),
                "mean": round(float(np.mean(arr)), 4),
                "std": round(float(np.std(arr, ddof=1)), 4) if len(arr) > 1 else 0.0,
                "min": round(float(np.min(arr)), 4),
                "max": round(float(np.max(arr)), 4),
            }
        stats_report[tool_name] = tool_stats

    return stats_report


def format_stats_table(stats: dict) -> str:
    """Format descriptive statistics as a readable table."""
    metrics = ["cai", "gc_content", "restriction_site_count", "gt_count",
               "constraints_satisfied", "constraint_violations", "execution_time_s"]
    metric_labels = {
        "cai": "CAI",
        "gc_content": "GC Content",
        "restriction_site_count": "Restriction Sites",
        "gt_count": "GT Dinucleotides",
        "constraints_satisfied": "Constraints Satisfied",
        "constraint_violations": "Constraint Violations",
        "execution_time_s": "Runtime (s)",
    }

    lines = []
    lines.append("=" * 90)
    lines.append("  Statistical Summary (12 genes × 4 methods)")
    lines.append("=" * 90)

    for metric in metrics:
        label = metric_labels.get(metric, metric)
        lines.append("")
        lines.append(f"  {label}")
        lines.append(f"  {'Method':<25s} {'Mean':>8s} {'Std':>8s} {'Min':>8s} {'Max':>8s} {'N':>4s}")
        lines.append("  " + "-" * 60)
        for tool_name in ["BC_constraint_first", "BC_cai_first", "SimpleCAI", "Random"]:
            if tool_name not in stats:
                continue
            m = stats[tool_name].get(metric, {})
            if m.get("n", 0) == 0:
                lines.append(f"  {tool_name:<25s}   (no data)")
                continue
            lines.append(
                f"  {tool_name:<25s} {m['mean']:>8.4f} {m['std']:>8.4f} "
                f"{m['min']:>8.4f} {m['max']:>8.4f} {m['n']:>4d}"
            )
    lines.append("")
    return "\n".join(lines)


# ============================================================================
# Part 3: Ablation Study
# ============================================================================

def run_ablation_study(
    genes: dict | None = None,
    enzymes: list[str] | None = None,
    gc_lo: float = 0.30,
    gc_hi: float = 0.70,
) -> list[dict]:
    """Run BioCompiler with each phase disabled to measure contribution."""
    genes = genes or GENE_PANEL
    enzymes = enzymes or DEFAULT_ENZYMES

    ablation_configs = [
        {"name": "Full_pipeline", "skip_phases": []},
        {"name": "Skip_Phase2_RS_removal", "skip_phases": [2]},
        {"name": "Skip_Phase3_cross_codon", "skip_phases": [3]},
        {"name": "Skip_Phase5_CpG_avoidance", "skip_phases": [5]},
        {"name": "Phase1_only", "skip_phases": [2, 3, 4, 5, 6, 7]},
    ]

    all_results = []
    for gene_name, (protein, organism) in genes.items():
        species = ORGANISM_TO_SPECIES.get(organism, "human")
        print(f"  Ablation: {gene_name} ({len(protein)} aa)...")

        for config in ablation_configs:
            result = _run_ablation_config(
                protein, organism, species, enzymes, gc_lo, gc_hi, config
            )
            result["gene"] = gene_name
            result["protein_length"] = len(protein)
            all_results.append(result)

    return all_results


def _run_ablation_config(
    protein: str,
    organism: str,
    species: str,
    enzymes: list[str],
    gc_lo: float,
    gc_hi: float,
    config: dict,
) -> dict:
    """Run BioCompiler with specific phases disabled."""
    skip = set(config["skip_phases"])
    t0 = time.perf_counter()
    try:
        opt = BioOptimizer(species=species, enzymes=enzymes, avoid_gt=True,
                           strategy="constraint_first")
        initial_seq = _build_best_codon_sequence(protein, species)
        seq = initial_seq.upper().strip()

        # Reset internal state
        opt._unavoidable_gt_positions = set()
        opt._applied_mutagenesis = []
        opt._original_protein = opt._translate(seq)

        # Phase 0: Max-CAI back-translation (always run)
        seq = opt._phase0_max_cai_backtranslate(seq)

        # Phase 1: Priority constraint resolution (always run)
        seq = opt._phase1_priority_constraint_resolution(seq)

        # Phase 2: Restriction site removal
        if 2 not in skip:
            seq = opt._phase2_remove_restriction_sites(seq)

        # Phase 3: Cross-codon constraint resolution
        if 3 not in skip:
            from .mutagenesis import MutagenesisReport
            seq, mut_report = opt._phase3_cross_codon_constraints(seq)
            seq, mut_report_35 = opt._phase35_within_codon_gt(seq)
            mut_report.proposals.extend(mut_report_35.proposals)
        else:
            from .mutagenesis import MutagenesisReport
            mut_report = MutagenesisReport()

        # Phase 4: Mutagenesis fallback
        if 4 not in skip:
            seq = opt._phase4_mutagenesis_fallback(seq, mut_report)

        # Phase 5: CpG island avoidance
        if 5 not in skip:
            seq = opt._phase5_avoid_cpg_islands(seq)

        # Phase 6: CAI hill climbing
        if 6 not in skip:
            seq = opt._phase6_cai_hill_climb(seq)

        # Phase 7: Re-optimization pass
        if 7 not in skip:
            seq = opt._phase7_reoptimize(seq)

        elapsed = time.perf_counter() - t0
        metrics = _compute_metrics(seq, protein, organism, species, enzymes, gc_lo, gc_hi)
        return {
            "ablation_config": config["name"],
            "sequence": seq,
            "execution_time_s": round(elapsed, 4),
            "success": True,
            **metrics,
        }
    except Exception as exc:
        elapsed = time.perf_counter() - t0
        logger.error("Ablation %s failed: %s", config["name"], exc)
        return {
            "ablation_config": config["name"],
            "sequence": "",
            "execution_time_s": round(elapsed, 4),
            "success": False,
            "error": str(exc),
            "cai": 0.0, "gc_content": 0.0,
            "restriction_site_count": -1, "gt_count": -1,
            "cpg_ratio": 0.0,
            "constraints_satisfied": 0, "constraint_violations": MAX_CONSTRAINTS,
            "max_constraints": MAX_CONSTRAINTS,
        }


def format_ablation_table(ablation_results: list[dict]) -> str:
    """Format ablation study results showing CAI impact of each phase."""
    valid = [r for r in ablation_results if r.get("success")]
    if not valid:
        return "No valid ablation results."

    configs = sorted(set(r["ablation_config"] for r in valid))
    lines = []
    lines.append("=" * 70)
    lines.append("  Ablation Study: CAI Impact of Each Phase")
    lines.append("=" * 70)
    lines.append(f"  {'Config':<30s} {'Mean CAI':>10s} {'Mean Violations':>16s} {'Mean Time (s)':>14s}")
    lines.append("  " + "-" * 70)

    # Compute full pipeline baseline for delta calculation
    full_cai = [r["cai"] for r in valid if r["ablation_config"] == "Full_pipeline"]
    full_cai_mean = np.mean(full_cai) if full_cai else 0

    for config in configs:
        vals = [r for r in valid if r["ablation_config"] == config]
        cai_mean = np.mean([r["cai"] for r in vals])
        viol_mean = np.mean([r["constraint_violations"] for r in vals])
        time_mean = np.mean([r["execution_time_s"] for r in vals])
        delta = cai_mean - full_cai_mean if config != "Full_pipeline" else 0
        delta_str = f" ({delta:+.4f})" if config != "Full_pipeline" else " (baseline)"
        lines.append(
            f"  {config:<30s} {cai_mean:>10.4f}{delta_str:<16s} "
            f"{viol_mean:>8.2f}        {time_mean:>8.4f}"
        )

    lines.append("")
    return "\n".join(lines)


# ============================================================================
# Part 4: Pareto Frontier Plot
# ============================================================================

def plot_pareto_frontier(
    results: list[dict],
    output_dir: Path,
) -> str:
    """Generate Pareto frontier plot: CAI (x) vs constraint_violations (y)."""
    output_dir.mkdir(parents=True, exist_ok=True)

    tool_styles = {
        "BC_constraint_first": {"color": "#e74c3c", "marker": "*", "s": 200, "zorder": 5, "label": "BC constraint_first"},
        "BC_cai_first": {"color": "#e67e22", "marker": "D", "s": 120, "zorder": 5, "label": "BC cai_first"},
        "SimpleCAI": {"color": "#2ecc71", "marker": "^", "s": 100, "zorder": 3, "label": "SimpleCAI"},
        "Random": {"color": "#95a5a6", "marker": "o", "s": 70, "zorder": 2, "label": "Random"},
    }

    # Collect per-gene data
    gene_data: dict[str, list] = {}
    for r in results:
        gene = r["gene"]
        if gene not in gene_data:
            gene_data[gene] = []
        if r.get("success"):
            gene_data[gene].append(r)

    # ---- Per-gene subplot grid ----
    n_genes = len(gene_data)
    n_cols = 4
    n_rows = math.ceil(n_genes / n_cols)

    fig, axes = plt.subplots(n_rows, n_cols, figsize=(5 * n_cols, 4 * n_rows))
    if n_genes == 1:
        axes = np.array([[axes]])
    axes_flat = axes.flatten()

    for idx, (gene, tools) in enumerate(gene_data.items()):
        if idx >= len(axes_flat):
            break
        ax = axes_flat[idx]
        seen = set()
        for t in tools:
            style = tool_styles.get(t["tool"],
                {"color": "gray", "marker": "x", "s": 40, "zorder": 1, "label": t["tool"]})
            lbl = style.get("label", t["tool"]) if t["tool"] not in seen else None
            ax.scatter(
                t["cai"], t["constraint_violations"],
                c=style["color"], marker=style["marker"], s=style["s"],
                label=lbl, zorder=style["zorder"],
                edgecolors="black", linewidths=0.5,
            )
            seen.add(t["tool"])

        ax.set_xlabel("CAI", fontsize=9)
        ax.set_ylabel("Constraint Violations", fontsize=9)
        ax.set_title(gene, fontsize=10, fontweight="bold")
        ax.set_xlim(0, 1.05)
        ax.set_ylim(-0.3, MAX_CONSTRAINTS + 0.5)
        ax.set_yticks(range(MAX_CONSTRAINTS + 1))
        ax.grid(True, alpha=0.3)

    # Hide empty subplots
    for idx in range(len(gene_data), len(axes_flat)):
        axes_flat[idx].set_visible(False)

    # Common legend
    handles, labels = [], []
    for tool_name, style in tool_styles.items():
        h = plt.scatter([], [], c=style["color"], marker=style["marker"],
                        s=style["s"] // 2, edgecolors="black", linewidths=0.5)
        handles.append(h)
        labels.append(style.get("label", tool_name))
    fig.legend(handles, labels, loc='upper center', ncol=4, fontsize=9,
               bbox_to_anchor=(0.5, 1.02))

    plt.tight_layout(rect=[0, 0, 1, 0.96])
    fpath = output_dir / "pareto_frontier.png"
    fig.savefig(fpath, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"  Saved: {fpath}")
    return str(fpath)


# ============================================================================
# Part 5: Ablation Plot
# ============================================================================

def plot_ablation(
    ablation_results: list[dict],
    output_dir: Path,
) -> str:
    """Generate ablation study bar chart."""
    output_dir.mkdir(parents=True, exist_ok=True)
    valid = [r for r in ablation_results if r.get("success")]
    if not valid:
        return ""

    configs = sorted(set(r["ablation_config"] for r in valid))
    metrics = ["cai", "constraint_violations", "execution_time_s"]
    metric_labels = ["CAI", "Constraint Violations", "Runtime (s)"]

    fig, axes = plt.subplots(1, 3, figsize=(18, 6))

    for ax_idx, (metric, label) in enumerate(zip(metrics, metric_labels)):
        ax = axes[ax_idx]
        means = []
        stds = []
        for config in configs:
            vals = [r[metric] for r in valid if r["ablation_config"] == config]
            means.append(np.mean(vals) if vals else 0)
            stds.append(np.std(vals, ddof=1) if len(vals) > 1 else 0)

        x = np.arange(len(configs))
        colors = plt.cm.Set2(np.linspace(0, 1, len(configs)))
        bars = ax.bar(x, means, yerr=stds, capsize=3, color=colors,
                      edgecolor="black", linewidth=0.5)
        ax.set_ylabel(label, fontsize=11)
        ax.set_xticks(x)
        short = [c.replace("Skip_Phase2_RS_removal", "-RS")
                 .replace("Skip_Phase3_cross_codon", "-CrossCodon")
                 .replace("Skip_Phase5_CpG_avoidance", "-CpG")
                 .replace("Full_pipeline", "Full")
                 .replace("Phase1_only", "Phase1")
                 for c in configs]
        ax.set_xticklabels(short, rotation=45, ha="right", fontsize=9)
        ax.set_title(f"Ablation: {label}", fontsize=12, fontweight="bold")
        ax.grid(True, alpha=0.3, axis="y")

        for bar, mean_val in zip(bars, means):
            ax.text(bar.get_x() + bar.get_width() / 2., bar.get_height() + 0.01,
                    f'{mean_val:.3f}', ha='center', va='bottom', fontsize=7)

    plt.tight_layout()
    fpath = output_dir / "ablation_study.png"
    fig.savefig(fpath, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"  Saved: {fpath}")
    return str(fpath)


# ============================================================================
# Output Persistence
# ============================================================================

def save_json(data: dict, filepath: Path) -> None:
    """Save data as JSON."""
    filepath.parent.mkdir(parents=True, exist_ok=True)
    with open(filepath, "w", encoding="utf-8") as f:
        json.dump(data, f, indent=2, default=str, ensure_ascii=False)
    print(f"  Saved: {filepath}")


def save_csv(results: list[dict], filepath: Path) -> None:
    """Save results as CSV (tabular summary)."""
    filepath.parent.mkdir(parents=True, exist_ok=True)
    columns = [
        "gene", "tool", "protein_length", "organism",
        "cai", "gc_content", "restriction_site_count", "gt_count",
        "cpg_ratio", "constraints_satisfied", "constraint_violations",
        "execution_time_s", "success",
    ]
    # Also include ablation_config if present
    has_ablation = any("ablation_config" in r for r in results)
    if has_ablation:
        columns.insert(2, "ablation_config")

    with open(filepath, "w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=columns, extrasaction="ignore")
        writer.writeheader()
        for r in results:
            row = {k: r.get(k, "") for k in columns}
            writer.writerow(row)
    print(f"  Saved: {filepath}")


def save_summary(text: str, filepath: Path) -> None:
    """Save text summary."""
    filepath.parent.mkdir(parents=True, exist_ok=True)
    with open(filepath, "w", encoding="utf-8") as f:
        f.write(text)
    print(f"  Saved: {filepath}")


# ============================================================================
# Main Entry Point
# ============================================================================

def main():
    """Run the comprehensive benchmark suite."""
    print("=" * 70)
    print("  BioCompiler Comprehensive Benchmark Suite v2.0")
    print("  12 genes × 4 methods + ablation + statistical analysis")
    print("=" * 70)
    print()

    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

    # ── Part 1: Multi-gene comparison ──────────────────────────
    print("Part 1: Multi-gene comparison (12 genes × 4 methods)...")
    benchmark_results = run_multi_gene_comparison()

    # ── Part 2: Statistical analysis ───────────────────────────
    print("\nPart 2: Statistical analysis...")
    stats = compute_statistics(benchmark_results)
    stats_table = format_stats_table(stats)
    print(stats_table)

    # ── Part 3: Ablation study ─────────────────────────────────
    print("\nPart 3: Ablation study...")
    ablation_results = run_ablation_study()
    ablation_table = format_ablation_table(ablation_results)
    print(ablation_table)

    # ── Part 4: Pareto frontier plot ───────────────────────────
    print("\nPart 4: Generating Pareto frontier plot...")
    pareto_file = plot_pareto_frontier(benchmark_results, OUTPUT_DIR)

    # ── Part 5: Ablation bar chart ─────────────────────────────
    print("\nPart 5: Generating ablation bar chart...")
    ablation_plot_file = plot_ablation(ablation_results, OUTPUT_DIR)

    # ── Save output files ──────────────────────────────────────
    print("\nSaving output files...")

    # comprehensive_results.json
    save_json({
        "timestamp": datetime.now(timezone.utc).isoformat(),
        "gene_panel": {k: {"protein_length": len(v[0]), "organism": v[1]}
                       for k, v in GENE_PANEL.items()},
        "benchmark_results": benchmark_results,
        "ablation_results": ablation_results,
        "statistics": stats,
    }, OUTPUT_DIR / "comprehensive_results.json")

    # comprehensive_results.csv (benchmark + ablation combined)
    all_rows = benchmark_results + ablation_results
    save_csv(all_rows, OUTPUT_DIR / "comprehensive_results.csv")

    # summary.txt
    summary = "\n".join([
        "=" * 70,
        "  BioCompiler Comprehensive Benchmark Summary",
        f"  Generated: {datetime.now(timezone.utc).isoformat()}",
        "=" * 70,
        "",
        f"Genes tested: {len(GENE_PANEL)}",
        f"Methods compared: BC_constraint_first, BC_cai_first, SimpleCAI, Random",
        f"Enzymes checked: {', '.join(DEFAULT_ENZYMES)}",
        "",
        stats_table,
        "",
        ablation_table,
        "",
        "Output files:",
        f"  - comprehensive_results.json",
        f"  - comprehensive_results.csv",
        f"  - pareto_frontier.png",
        f"  - ablation_study.png",
        f"  - summary.txt",
        "",
    ])
    save_summary(summary, OUTPUT_DIR / "summary.txt")

    print("\n" + "=" * 70)
    print("  Benchmark complete!")
    print(f"  Results saved to: {OUTPUT_DIR}")
    print("=" * 70)


if __name__ == "__main__":
    # Ensure src/ is prioritized over any root-level biocompiler/ package
    import sys
    src_dir = str(Path(__file__).resolve().parent.parent)
    if src_dir not in sys.path:
        sys.path.insert(0, src_dir)
    logging.basicConfig(level=logging.INFO, format="%(levelname)s: %(message)s")
    main()
