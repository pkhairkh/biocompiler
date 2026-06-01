"""
BioCompiler Benchmark Module v7.0.0
=====================================
Built-in benchmark sequences and performance measurement for BioCompiler.

Provides:
  - run_benchmark(): Optimize eGFP, mCherry, LacZ with human/ecoli and print a table
  - compare_tools(): Theoretical feature comparison table vs other tools
"""

import time
from typing import Dict, List, Tuple

try:
    from .optimizer import BioOptimizer
except ImportError:
    from .optimization import BioOptimizer
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

# eGFP (Enhanced Green Fluorescent Protein) — 717 bp coding sequence
# Source: pEGFP-N1 (Clontech), GenBank accession U55763 region
EGFP_DNA = (
    "ATGGTGAGCAAGGGCGAGGAGCTGTTCACCGGGGTGGTGCCCATCCTGGTCGAGCTGGACGGCGACGTAAACGGCCACAAGTTCAGCGTGTCCGGCGAGGGCGAGGGCGATGCCACCTACGGCAAGCTGACCCTGAAGTTCATCTGCACCACCGGCAAGCTGCCCGTGCCCTGGCCCACCCTCGTGACCACCCTGACCTACGGCGTGCAGTGCTTCAGCCGCTACCCCGACCACATGAAGCAGCACGACTTCTTCAAGTCCGCCATGCCCGAAGGCTACGTCCAGGAGCGCACCATCTTCTTCAAGGACGACGGCAACTACAAGACCCGCGCCGAGGTGAAGTTCGAGGGCGACACCCTGGTGAACCGCATCGAGCTGAAGGGCATCGACTTCAAGGAGGACGGCAACATCCTGGGGCACAAGCTGGAGTACAACTACAACAGCCACAACGTCTATATCATGGCCGACAAGCAGAAGAACGGCATCAAGGTGAACTTCAAGATCCGCCACAACATCGAGGACGGCAGCGTGCAGCTCGCCGACCACTACCAGCAGAACACCCCCATCGGCGACGGCCCCGTGCTGCTGCCCGACAACCACTACCTGAGCACCCAGTCCGCCCTGAGCAAAGACCCCAACGAGAAGCGCGATCACATGGTCCTGCTGGAGTTCGTGACCGCCGCCGGGATCACTCTCGGCATGGACGAGCTGTACAAGTAA"
)

# mCherry (Red Fluorescent Protein) — 711 bp coding sequence
# Source: pmCherry (Clontech), derived from mRFP1
MCHERRY_DNA = (
    "ATGGTGAGCAAGGGCGAGGAGGATAACATGGCCATCATCAAGGAGTTCATGCGCTTCAAGGTGCACATGGAGGGCTCCGTGAACGGCCACGAGTTCGAGATCGAGGGCGAGGGCGAGGGCCGCCCCTACGAGGGCACCCAGACCGCCAAGCTGAAGGTGACCAAGGGTGGCCCCCTGCCCTTCGCCTGGGACATCCTGTCCCCTCAGTTCATGTACGGCTCCAAGGCCTACGTGAAGCACCCCGCCGACATCCCCGACTACTTGAAGCTGTCCTTCCCCGAGGGCTTCAAGTGGGAGCGCGTGATGAACTTCGAGGACGGCGGCGTGGTGACCGTGACCCAGGACTCCTCCCTGCAGGACGGCGAGTTCATCTACAAGGTGAAGCTGCGCGGCACCAACTTCCCCTCCGACGGCCCCGTAATGCAGAAGAAGACCATGGGCTGGGAGGCCTCCTCCGAGCGGATGTACCCCGAGGACGGCGCCCTGAAGGGCGAGATCAAGCAGCGGCTGAAGCTGAAGGACGGCGGCCACTACGACGCTGAGGTCAAGACCACCTACAAGGCCAAGAAGCCCGTGCAGCTGCCCGGCGCCTACAACGTCAACATCAAGTTGGACATCACCTCCCACAACGAGGACTACACCATCGTGGAACAGTACGAACGCGCCGAGGGCCGCCACTCCACCGGCGGCATGGACGAGCTGTACAAGTAA"
)

# LacZ (beta-galactosidase, N-terminal 720 bp fragment) — 720 bp
# Source: E. coli lacZ gene, first 240 codons
LACZ_DNA = (
    "ATGACCATGATTACGCCAAGCTATTTAGGTGACACTATAGAATACTCAAGCTATGCATCCAACGCGTTGGGAGCTCTCCCATATGGTCGACCTGCAGGCGGCCGCACTAGTGATTACGCCAAGCTATTTAGGTGACACTATAGAATACTCAAGCTATGCATCCAACGCGTTGGGAGCTCTCCCATATGGTCGACCTGCAGGCGGCCGCACTAGTGATTATGCCTGCAGGTCGACTCTAGAGGATCCCGGGTACCGAGCTCGAATTCACTGGCCGTCGTTTTACAACGTCGTGACTGGGAAAACCCTGGCGTTACCCAACTTAATCGCCTTGCAGCACATCCCCCTTTCGCCAGCTGGCGTAATAGCGAAGAGGCCCGCACCGATCGCCCTTCCCAACAGTTGCGCAGCCTGAATGGCGAATGGCGCTTTGCCTGGTTTCCGGCACCAGAAGCGGTGCCGGAAAGCTGGCTGGAGTGCGATCTTCCTGAGGCCGATACTGTCGTCGTCCCCTCAAACTGGCAGATGCACGGTTACGATGCGCCCATCTACACCAACGTGACCTATCCCATTACGGTCAATCCGCCGTTTGTTCCCACGGAGAATCCGACGGGTTGTTACTCGCTCACATTTAATGTTGATGAAAGCTGGCTACAGGAAGGCCAGACGCGAATTATTTTTGATGGCGTTAACTCGGCGTTTCATCTGTGGTGCAACGGGCGCTGGGTCGGTTACGGCCAGGACAGTCGTTTGCCGTCTGAATTTGACCTGAGCGCATTTTTACGCGCCGGAGAAAACCGCCTCGCGGTGATGGTGCTGCGTTGGAGTGACGGCAGTTATCTGGAAGATCAGGATATGTGGCGGATGAGCGGCATTTTCCGTGACGTCTCGTTGCTGCATAAACCGACTACACAAATCAGCGATTTCCATGTTGCCACTCGCTTTAATGATGATTTCAGCCGCGCTGTACTGGAGGCTGAAGTTCAGATGTGCGGCGAGTTGCGTGACTACCTACGGGTAACAGTTTCTTTATGGCAGGGTGAAACGCAGGTCGCCAGCGGCACCGCGCCTTTCGGCGGTGAAATTATCGATGAGCGTGGTGGTTATGCCGATCGC"
)


def _compute_cai(seq: str, species_cai: Dict[str, float]) -> float:
    """Compute the geometric mean CAI for a sequence."""
    import math
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
) -> None:
    """Run built-in benchmarks for eGFP, mCherry, and LacZ.

    Optimizes each gene with both human and ecoli species tables,
    then prints a formatted table of results.

    Args:
        enzymes: List of restriction enzymes to avoid
        splice_low: Low splice threshold
        splice_high: High splice threshold
    """
    if enzymes is None:
        enzymes = ["EcoRI", "BamHI", "HindIII", "XhoI"]

    genes = {
        "eGFP": EGFP_DNA,
        "mCherry": MCHERRY_DNA,
        "LacZ": LACZ_DNA,
    }

    species_list = ["human", "ecoli"]

    print()
    print("=" * 100)
    print("  BioCompiler v7.0.0 — Built-in Benchmark")
    print("=" * 100)
    print(f"  Enzymes avoided: {', '.join(enzymes)}")
    print(f"  Splice thresholds: low={splice_low}, high={splice_high}")
    print("=" * 100)

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
    print("=" * 100)
    print()


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
            },
        },
    ]

    feature_names = list(tools[0]["features"].keys())

    print()
    print("=" * 80)
    print("  BioCompiler v7.0.0 — Tool Comparison (Theoretical)")
    print("=" * 80)
    print()

    # Header
    name_col = 20
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
    print("    - BioCompiler:   Full pipeline with formal GOLD/SILVER/BRONZE certificates")
    print("=" * 80)
    print()


# Alias for compatibility with tests that use run_benchmarks (plural)
run_benchmarks = run_benchmark
