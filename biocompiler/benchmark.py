"""
BioCompiler Benchmark Module v8.0.0
=====================================
Built-in benchmark sequences and performance measurement for BioCompiler.

Provides:
  - run_benchmark(): Optimize eGFP, mCherry, LacZ, Insulin, HBB across species
  - run_extended_benchmark(): Full extended benchmark with CAI-only comparison
  - compare_tools(): Theoretical feature comparison table vs other tools
"""

import math
import time
from typing import Dict, List, Tuple

from .optimizer import BioOptimizer
from .type_system import (
    CODON_TABLE,
    CertLevel,
    PredicateResult,
    check_no_avoidable_gt,
    check_no_cpg_island,
    check_no_restriction_site,
)
from .certificates import compute_certificate
from .species import SPECIES

# ────────────────────────────────────────────────────────────
# Built-in gene sequences (standard reference sequences)
# ────────────────────────────────────────────────────────────

# eGFP (Enhanced Green Fluorescent Protein) — 720 bp coding sequence
# Source: pEGFP-N1 (Clontech), GenBank accession U55763 region
EGFP_DNA = (
    "ATGGTGAGCAAGGGCGAGGAGCTGTTCACCGGGGTGGTGCCCATCCTGGTCGAGCTGGACGGCGACGTAAACGGCCACAAGTTCAGCGTGTCCGGCGAGGGCGAGGGCGATGCCACCTACGGCAAGCTGACCCTGAAGTTCATCTGCACCACCGGCAAGCTGCCCGTGCCCTGGCCCACCCTCGTGACCACCCTGACCTACGGCGTGCAGTGCTTCAGCCGCTACCCCGACCACATGAAGCAGCACGACTTCTTCAAGTCCGCCATGCCCGAAGGCTACGTCCAGGAGCGCACCATCTTCTTCAAGGACGACGGCAACTACAAGACCCGCGCCGAGGTGAAGTTCGAGGGCGACACCCTGGTGAACCGCATCGAGCTGAAGGGCATCGACTTCAAGGAGGACGGCAACATCCTGGGGCACAAGCTGGAGTACAACTACAACAGCCACAACGTCTATATCATGGCCGACAAGCAGAAGAACGGCATCAAGGTGAACTTCAAGATCCGCCACAACATCGAGGACGGCAGCGTGCAGCTCGCCGACCACTACCAGCAGAACACCCCCATCGGCGACGGCCCCGTGCTGCTGCCCGACAACCACTACCTGAGCACCCAGTCCGCCCTGAGCAAAGACCCCAACGAGAAGCGCGATCACATGGTCCTGCTGGAGTTCGTGACCGCCGCCGGGATCACTCTCGGCATGGACGAGCTGTACAAGTAA"
)

# mCherry (Red Fluorescent Protein) — 711 bp coding sequence
# Source: pmCherry (Clontech), derived from mRFP1
MCHERRY_DNA = (
    "ATGGTGAGCAAGGGCGAGGAGGATAACATGGCCATCATCAAGGAGTTCATGCGCTTCAAGGTGCACATGGAGGGCTCCGTGAACGGCCACGAGTTCGAGATCGAGGGCGAGGGCGAGGGCCGCCCCTACGAGGGCACCCAGACCGCCAAGCTGAAGGTGACCAAGGGTGGCCCCCTGCCCTTCGCCTGGGACATCCTGTCCCCTCAGTTCATGTACGGCTCCAAGGCCTACGTGAAGCACCCCGCCGACATCCCCGACTACTTGAAGCTGTCCTTCCCCGAGGGCTTCAAGTGGGAGCGCGTGATGAACTTCGAGGACGGCGGCGTGGTGACCGTGACCCAGGACTCCTCCCTGCAGGACGGCGAGTTCATCTACAAGGTGAAGCTGCGCGGCACCAACTTCCCCTCCGACGGCCCCGTAATGCAGAAGAAGACCATGGGCTGGGAGGCCTCCTCCGAGCGGATGTACCCCGAGGACGGCGCCCTGAAGGGCGAGATCAAGCAGCGGCTGAAGCTGAAGGACGGCGGCCACTACGACGCTGAGGTCAAGACCACCTACAAGGCCAAGAAGCCCGTGCAGCTGCCCGGCGCCTACAACGTCAACATCAAGTTGGACATCACCTCCCACAACGAGGACTACACCATCGTGGAACAGTACGAACGCGCCGAGGGCCGCCACTCCACCGGCGGCATGGACGAGCTGTACAAGTAA"
)

# LacZ (beta-galactosidase, N-terminal 867 bp fragment) — 867 bp
# Source: E. coli lacZ gene, GenBank V00296, first 289 codons
LACZ_DNA = (
    "ATGACCATGATTACGGATTCACTGGCCGTCGTTTTACAACGTCGTGACTGGGAAAACCCTGGCGTTACCCAACTTAATCGCCTTGCAGCACATCCCCCTTTCGCCAGCTGGCGTAATAGCGAAGAGGCCCGCACCGATCGCCCTTCCCAACAGTTGCGCAGCCTGAATGGCGAATGGCGCTTTGCCTGGTTTCCGGCACCAGAAGCGGTGCCGGAAAGCTGGCTGGAGTGCGATCTTCCTGAGGCCGATACTGTCGTCGTCCCCTCAAACTGGCAGATGCACGGTTACGATGCGCCCATCTACACCAACGTGACCTATCCCATTACGGTCAATCCGCCGTTTGTTCCCACGGAGAATCCGACGGGTTGTTACTCGCTCACATTTAATGTTGATGAAAGCTGGCTACAGGAAGGCCAGACGCGAATTATTTTTGATGGCGTTAACTCGGCGTTTCATCTGTGGTGCAACGGGCGCTGGGTCGGTTACGGCCAGGACAGTCGTTTGCCGTCTGAATTTGACCTGAGCGCATTTTTACGCGCCGGAGAAAACCGCCTCGCGGTGATGGTGCTGCGTTGGAGTGACGGCAGTTATCTGGAAGATCAGGATATGTGGCGGATGAGCGGCATTTTCCGTGACGTCTCGTTGCTGCATAAACCGACTACACAAATCAGCGATTTCCATGTTGCCACTCGCTTTAATGATGATTTCAGCCGCGCTGTACTGGAGGCTGAAGTTCAGATGTGCGGCGAGTTGCGTGACTACCTACGGGTAACAGTTTCTTTATGGCAGGGTGAAACGCAGGTCGCCAGCGGCACCGCGCCTTTCGGCGGTGAAATTATCGATGAGCGTGGTGGTTATGCCGATCGC"
)

# Human Insulin (preproinsulin) — 333 bp coding sequence
# Source: Human INS gene, UniProt P01308, back-translated from 110 aa protein
INSULIN_DNA = (
    "ATGGCCCTGTGGATGAGGCTGCTGCCCCTGCTGGCCCTGCTGGCCCTGTGGGGCCCCGACCCCGCCGCCGCCTTCGTGAACCAGCACCTGTGCGGCAGCCACCTGGTGGAGGCCCTGTACCTGGTGTGCGGCGAGAGGGGCTTCTTCTACACCCCCAAGACCAGGAGGGAGGCCGAGGACCTGCAGGTGGGCCAGGTGGAGCTGGGCGGCGGCCCCGGCGCCGGCAGCCTGCAGCCCCTGGCCCTGGAGGGCAGCCTGCAGAAGAGGGGCATCGTGGAGCAGTGCTGCACCAGCATCTGCAGCCTGTACCAGCTGGAGAACTACTGCAACTAA"
)

# Human beta-globin (HBB) — 444 bp coding sequence
# Source: Human HBB gene, GenBank NM_000518.5
HBB_DNA = (
    "ATGGTGCACCTGACTCCTGAGGAGAAGTCTGCCGTTACTGCCCTGTGGGGCAAGGTGAACGTGGATGAAGTTGGTGGTGAGGCCCTGGGCAGGCTGCTGGTGGTCTACCCTTGGACCCAGAGGTTCTTTGAGTCCTTTGGGGATCTGTCCACTCCTGATGCTGTTATGGGCAACCCTAAGGTGAAGGCTCATGGCAAGAAAGTGCTCGGTGCCTTTAGTGATGGCCTGGCTCACCTGGACAACCTCAAGGGCACCTTTGCCACACTGAGTGAGCTGCACTGTGACAAGCTGCACGTGGATCCTGAGAACTTCAGGCTCCTGGGCAACGTGCTGGTCTGTGTGCTGGCCCATCACTTTGGCAAAGAATTCACCCCACCAGTGCAGGCTGCCTATCAGAAAGTGGTGGCTGGTGTGGCTAATGCCCTGGCCCACAAGTATCACTAA"
)


def _compute_cai(seq: str, species_cai: Dict[str, float]) -> float:
    """Compute the geometric mean CAI for a sequence."""
    if not seq or len(seq) < 3:
        return 0.0
    log_sum = 0.0
    count = 0
    for i in range(0, len(seq) - 2, 3):
        codon = seq[i:i+3]
        cai = species_cai.get(codon, 0.0)
        if cai <= 0:
            cai = 0.001  # avoid log(0)
        log_sum += math.log(cai)
        count += 1
    if count == 0:
        return 0.0
    return math.exp(log_sum / count)


def _count_gt(seq: str) -> int:
    """Count GT dinucleotides in a sequence."""
    return sum(1 for i in range(len(seq) - 1) if seq[i:i+2] == "GT")


def _count_cpg_ratio(seq: str) -> float:
    """Compute CpG Obs/Exp ratio for the full sequence."""
    c = seq.count("C")
    g = seq.count("G")
    cg = sum(1 for i in range(len(seq) - 1) if seq[i:i+2] == "CG")
    expected = (c * g) / len(seq) if len(seq) > 0 else 0
    return cg / expected if expected > 0 else 0.0


def _count_restriction_sites(seq: str, enzymes: List[str]) -> int:
    """Count total restriction sites in sequence."""
    from .restriction_sites import get_recognition_site
    total = 0
    for enzyme in enzymes:
        site = get_recognition_site(enzyme)
        if site is None:
            continue
        pos = seq.find(site)
        while pos != -1:
            total += 1
            pos = seq.find(site, pos + 1)
    return total


def run_benchmark(
    enzymes: List[str] = None,
    splice_low: float = 3.0,
    splice_high: float = 6.0,
    species_list: List[str] = None,
    genes: Dict[str, str] = None,
) -> List[Dict]:
    """Run built-in benchmarks for standard genes across species.

    Args:
        enzymes: List of restriction enzymes to avoid
        splice_low: Low splice threshold
        splice_high: High splice threshold
        species_list: Species to test (default: human, ecoli)
        genes: Dict of gene_name -> DNA sequence (default: eGFP, mCherry, LacZ)

    Returns:
        List of result dicts for further analysis
    """
    if enzymes is None:
        enzymes = ["EcoRI", "BamHI", "HindIII", "XhoI"]
    if species_list is None:
        species_list = ["human", "ecoli"]
    if genes is None:
        genes = {
            "eGFP": EGFP_DNA,
            "mCherry": MCHERRY_DNA,
            "LacZ": LACZ_DNA,
        }

    print()
    print("=" * 110)
    print("  BioCompiler v8.1.0 — Built-in Benchmark")
    print("=" * 110)
    print(f"  Genes: {', '.join(genes.keys())}")
    print(f"  Species: {', '.join(species_list)}")
    print(f"  Enzymes avoided: {', '.join(enzymes)}")
    print(f"  Splice thresholds: low={splice_low}, high={splice_high}")
    print("=" * 110)

    results: List[Dict] = []

    for gene_name, gene_seq in genes.items():
        for species in species_list:
            opt = BioOptimizer(
                species=species,
                enzymes=enzymes,
                splice_low=splice_low,
                splice_high=splice_high,
                avoid_gt=True,
            )

            # Measure optimization time
            t0 = time.perf_counter()
            optimized, pred_results, cert_text = opt.optimize(gene_seq)
            t1 = time.perf_counter()
            elapsed = t1 - t0

            species_cai = SPECIES.get(species, SPECIES["ecoli"])

            cai_before = _compute_cai(gene_seq, species_cai)
            cai_after = _compute_cai(optimized, species_cai)
            gt_before = _count_gt(gene_seq)
            gt_after = _count_gt(optimized)
            cpg_ratio = _count_cpg_ratio(optimized)
            rs_count = _count_restriction_sites(optimized, enzymes)
            cert_level = compute_certificate(pred_results)

            results.append({
                "gene": gene_name,
                "species": species,
                "length": len(optimized),
                "cai_before": cai_before,
                "cai_after": cai_after,
                "gt_before": gt_before,
                "gt_after": gt_after,
                "cpg_ratio": cpg_ratio,
                "rs_count": rs_count,
                "cert": cert_level.value,
                "time_ms": elapsed * 1000,
            })

    # Print formatted table
    print()
    header = (
        f"{'Gene':<10} {'Species':<8} {'Length':>6} "
        f"{'CAI before':>10} {'CAI after':>10} "
        f"{'GT before':>10} {'GT after':>9} "
        f"{'CpG ratio':>10} {'RS sites':>9} "
        f"{'Cert':>8} {'Time(ms)':>9}"
    )
    print(header)
    print("-" * len(header))

    for r in results:
        row = (
            f"{r['gene']:<10} {r['species']:<8} {r['length']:>6} "
            f"{r['cai_before']:>10.4f} {r['cai_after']:>10.4f} "
            f"{r['gt_before']:>10d} {r['gt_after']:>9d} "
            f"{r['cpg_ratio']:>10.3f} {r['rs_count']:>9d} "
            f"{r['cert']:>8} {r['time_ms']:>9.1f}"
        )
        print(row)

    print()
    print("  Legend: CAI = Codon Adaptation Index, GT = GT dinucleotides,")
    print("          CpG ratio = Obs/Exp CG, RS = Restriction sites, Cert = Certificate level")
    print("=" * 110)
    print()

    return results


def run_extended_benchmark(
    enzymes: List[str] = None,
    splice_low: float = 3.0,
    splice_high: float = 6.0,
) -> List[Dict]:
    """Run extended benchmark with all species, all genes, and CAI-only comparison.

    Tests 5 genes x 4 species with both GT-aware and CAI-only modes,
    then prints a comprehensive comparison table showing:
    - Full pipeline (GT-aware) results
    - CAI-only mode results (no GT avoidance, for comparison)
    - CAI retention: ratio of GT-aware CAI to CAI-only CAI
    """
    if enzymes is None:
        enzymes = ["EcoRI", "BamHI", "HindIII", "XhoI"]

    all_genes = {
        "eGFP": EGFP_DNA,
        "mCherry": MCHERRY_DNA,
        "LacZ": LACZ_DNA,
        "Insulin": INSULIN_DNA,
        "HBB": HBB_DNA,
    }

    all_species = ["human", "ecoli", "yeast", "cho"]

    print()
    print("=" * 120)
    print("  BioCompiler v8.1.0 — Extended Benchmark (5 Genes x 4 Species)")
    print("=" * 120)
    print(f"  Genes: {', '.join(all_genes.keys())}")
    print(f"  Species: {', '.join(all_species)}")
    print(f"  Enzymes avoided: {', '.join(enzymes)}")
    print(f"  Modes: GT-aware (full pipeline) + CAI-only (no GT avoidance)")
    print("=" * 120)

    results: List[Dict] = []

    for gene_name, gene_seq in all_genes.items():
        for species in all_species:
            species_cai = SPECIES.get(species, SPECIES["ecoli"])

            # Mode 1: Full pipeline (GT-aware)
            opt_gt = BioOptimizer(
                species=species,
                enzymes=enzymes,
                splice_low=splice_low,
                splice_high=splice_high,
                avoid_gt=True,
            )
            t0 = time.perf_counter()
            opt_seq_gt, pred_gt, cert_gt = opt_gt.optimize(gene_seq)
            t1 = time.perf_counter()

            cai_before = _compute_cai(gene_seq, species_cai)
            cai_gt = _compute_cai(opt_seq_gt, species_cai)
            gt_gt = _count_gt(opt_seq_gt)
            cpg_gt = _count_cpg_ratio(opt_seq_gt)
            rs_gt = _count_restriction_sites(opt_seq_gt, enzymes)
            cert_gt_level = compute_certificate(pred_gt)

            # Mode 2: CAI-only (no GT avoidance)
            opt_cai = BioOptimizer(
                species=species,
                enzymes=enzymes,
                splice_low=splice_low,
                splice_high=splice_high,
                avoid_gt=False,
            )
            t2 = time.perf_counter()
            opt_seq_cai, pred_cai, cert_cai = opt_cai.optimize(gene_seq)
            t3 = time.perf_counter()

            cai_cai = _compute_cai(opt_seq_cai, species_cai)
            gt_cai = _count_gt(opt_seq_cai)
            cpg_cai = _count_cpg_ratio(opt_seq_cai)
            rs_cai = _count_restriction_sites(opt_seq_cai, enzymes)

            # CAI retention: how much CAI does GT-aware mode retain vs CAI-only
            cai_retention = cai_gt / cai_cai if cai_cai > 0 else 0.0

            results.append({
                "gene": gene_name,
                "species": species,
                "length": len(opt_seq_gt),
                "cai_before": cai_before,
                "cai_gt": cai_gt,
                "cai_cai": cai_cai,
                "cai_retention": cai_retention,
                "gt_before": _count_gt(gene_seq),
                "gt_gt": gt_gt,
                "gt_cai": gt_cai,
                "cpg_gt": cpg_gt,
                "cpg_cai": cpg_cai,
                "rs_gt": rs_gt,
                "rs_cai": rs_cai,
                "cert": cert_gt_level.value,
                "time_gt_ms": (t1 - t0) * 1000,
                "time_cai_ms": (t3 - t2) * 1000,
            })

    # ─── Table 1: GT-aware full pipeline ─────────────────────
    print()
    print("  ┌─────────────────────────────────────────────────────────────────────────────────────────────────┐")
    print("  │ TABLE 1: Full Pipeline (GT-aware, 7-phase)                                                     │")
    print("  └─────────────────────────────────────────────────────────────────────────────────────────────────┘")
    print()
    header1 = (
        f"{'Gene':<10} {'Species':<8} {'Length':>6} "
        f"{'CAI bef':>8} {'CAI aft':>8} "
        f"{'GT aft':>7} {'CpG':>6} "
        f"{'RS':>3} {'Cert':>8} {'Time':>7}"
    )
    print(header1)
    print("-" * len(header1))
    for r in results:
        row = (
            f"{r['gene']:<10} {r['species']:<8} {r['length']:>6} "
            f"{r['cai_before']:>8.4f} {r['cai_gt']:>8.4f} "
            f"{r['gt_gt']:>7d} {r['cpg_gt']:>6.3f} "
            f"{r['rs_gt']:>3d} {r['cert']:>8} {r['time_gt_ms']:>6.0f}ms"
        )
        print(row)

    # ─── Table 2: CAI-only mode (no GT avoidance) ─────────────
    print()
    print("  ┌─────────────────────────────────────────────────────────────────────────────────────────────────┐")
    print("  │ TABLE 2: CAI-only Mode (no GT avoidance)                                                        │")
    print("  └─────────────────────────────────────────────────────────────────────────────────────────────────┘")
    print()
    header2 = (
        f"{'Gene':<10} {'Species':<8} "
        f"{'CAI aft':>8} {'GT aft':>7} {'CpG':>6} "
        f"{'RS':>3} {'Time':>7}"
    )
    print(header2)
    print("-" * len(header2))
    for r in results:
        row = (
            f"{r['gene']:<10} {r['species']:<8} "
            f"{r['cai_cai']:>8.4f} {r['gt_cai']:>7d} {r['cpg_cai']:>6.3f} "
            f"{r['rs_cai']:>3d} {r['time_cai_ms']:>6.0f}ms"
        )
        print(row)

    # ─── Table 3: CAI Retention Analysis ─────────────────────
    print()
    print("  ┌─────────────────────────────────────────────────────────────────────────────────────────────────┐")
    print("  │ TABLE 3: CAI Retention (GT-aware CAI / CAI-only CAI)                                            │")
    print("  │   Shows how much CAI the GT-aware pipeline retains compared to the CAI-only upper bound         │")
    print("  └─────────────────────────────────────────────────────────────────────────────────────────────────┘")
    print()
    header3 = (
        f"{'Gene':<10} {'Species':<8} "
        f"{'CAI GT':>8} {'CAI only':>8} {'Retain%':>8} {'GT saved':>9}"
    )
    print(header3)
    print("-" * len(header3))
    for r in results:
        gt_saved = r['gt_before'] - r['gt_gt']
        row = (
            f"{r['gene']:<10} {r['species']:<8} "
            f"{r['cai_gt']:>8.4f} {r['cai_cai']:>8.4f} "
            f"{r['cai_retention']*100:>7.1f}% {gt_saved:>9d}"
        )
        print(row)

    # ─── Summary statistics ──────────────────────────────────
    avg_retention = sum(r['cai_retention'] for r in results) / len(results) if results else 0
    avg_cai_gt = sum(r['cai_gt'] for r in results) / len(results) if results else 0
    avg_cai_cai = sum(r['cai_cai'] for r in results) / len(results) if results else 0
    total_gt_saved = sum(r['gt_before'] - r['gt_gt'] for r in results)
    total_gt_remaining = sum(r['gt_gt'] for r in results)
    gold_count = sum(1 for r in results if r['cert'] == 'GOLD')
    silver_count = sum(1 for r in results if r['cert'] == 'SILVER')
    bronze_count = sum(1 for r in results if r['cert'] == 'BRONZE')

    print()
    print("  ┌─────────────────────────────────────────────────────────────────────────────────────────────────┐")
    print("  │ SUMMARY                                                                                          │")
    print("  └─────────────────────────────────────────────────────────────────────────────────────────────────┘")
    print(f"    Average CAI retention:     {avg_retention*100:.1f}%")
    print(f"    Average CAI (GT-aware):    {avg_cai_gt:.4f}")
    print(f"    Average CAI (CAI-only):    {avg_cai_cai:.4f}")
    print(f"    Total GTs eliminated:      {total_gt_saved}")
    print(f"    Total GTs remaining:       {total_gt_remaining} (unavoidable)")
    print(f"    Certificates: {gold_count} GOLD, {silver_count} SILVER, {bronze_count} BRONZE")
    print("=" * 120)
    print()

    return results


def compare_tools() -> None:
    """Print a theoretical feature comparison table: BioCompiler vs other tools.

    Based on published capabilities:
      - GeneOptimizer (GeneArt/Thermo): CAI optimization, no GT/CpG avoidance
      - DNAworks: Restriction site avoidance, CAI, no GT/CpG
      - OPTIMIZER: CAI optimization only, no GT/CpG/restriction
    """
    tools = [
        {
            "name": "BioCompiler",
            "features": {
                "CAI optimization": True,
                "GT avoidance": True,
                "CpG avoidance": True,
                "Restriction sites": True,
                "Formal certification": True,
                "Cross-codon awareness": True,
                "Dual-threshold splice": True,
                "CAI-boost re-pass": True,
                "Multi-codon look-ahead": True,
            },
        },
        {
            "name": "GeneOptimizer",
            "features": {
                "CAI optimization": True,
                "GT avoidance": False,
                "CpG avoidance": False,
                "Restriction sites": True,
                "Formal certification": False,
                "Cross-codon awareness": False,
                "Dual-threshold splice": False,
                "CAI-boost re-pass": False,
                "Multi-codon look-ahead": False,
            },
        },
        {
            "name": "DNAworks",
            "features": {
                "CAI optimization": True,
                "GT avoidance": False,
                "CpG avoidance": False,
                "Restriction sites": True,
                "Formal certification": False,
                "Cross-codon awareness": False,
                "Dual-threshold splice": False,
                "CAI-boost re-pass": False,
                "Multi-codon look-ahead": False,
            },
        },
        {
            "name": "OPTIMIZER",
            "features": {
                "CAI optimization": True,
                "GT avoidance": False,
                "CpG avoidance": False,
                "Restriction sites": False,
                "Formal certification": False,
                "Cross-codon awareness": False,
                "Dual-threshold splice": False,
                "CAI-boost re-pass": False,
                "Multi-codon look-ahead": False,
            },
        },
    ]

    feature_names = list(tools[0]["features"].keys())

    print()
    print("=" * 90)
    print("  BioCompiler v8.1.0 — Tool Comparison (Theoretical)")
    print("=" * 90)
    print()

    # Header
    name_col = 22
    feat_col = 8
    header = f"{'Feature':<{name_col}}"
    for tool in tools:
        header += f" {tool['name']:^{feat_col}}"
    print(header)
    print("-" * len(header))

    # Rows
    for feat in feature_names:
        row = f"{feat:<{name_col}}"
        for tool in tools:
            val = tool["features"][feat]
            mark = "\u2713" if val else "\u2717"
            row += f" {mark:^{feat_col}}"
        print(row)

    print()
    print("  \u2713 = Supported   \u2717 = Not supported")
    print()
    print("  Notes:")
    print("    - GeneOptimizer: CAI + restriction sites, but no GT/CpG avoidance")
    print("    - DNAworks:      CAI + restriction sites, but no GT/CpG/formal cert")
    print("    - OPTIMIZER:     CAI only, no constraint handling or certification")
    print("    - BioCompiler:   7-phase pipeline with CAI-boost, formal GOLD/SILVER/BRONZE certificates")
    print("=" * 90)
    print()
