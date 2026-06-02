#!/usr/bin/env python3
"""
BioCompiler Real Dataset Analysis
==================================
Run BioCompiler v7.0.0 against real biological datasets with full
type-directed mutagenesis. This script produces a detailed report
of optimizer performance, predicate satisfaction, and mutagenesis
impact on proteins that are impossible at codon level.

This is a standalone analysis/reporting script that:
  1. Runs full dataset validation (all organisms, all genes)
  2. Runs the benchmark suite (HBB, INS, EGFP)
  3. Performs detailed per-gene optimization analysis
  4. Analyzes GT-mandatory amino acids (Valine, Cysteine)
  5. Deep-dives HBB as the key mutagenesis proof case
  6. Cross-organism optimization comparison
  7. Restriction site elimination analysis
  8. Mutagenesis impact assessment on all failing genes
  9. Summary statistics by organism

Usage:
    python scripts/run_real_dataset.py
"""

import sys
import os
import time
import json

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "src"))

from biocompiler import (
    optimize_sequence,
    compute_cai,
    gc_content,
    translate,
    evaluate_all_predicates,
    type_directed_mutagenesis,
    find_unrepairable_cryptic_donors,
    find_unrepairable_cryptic_acceptors,
    BLOSUM62,
    GT_MANDATORY_AAS,
    run_dataset_validation,
)
from biocompiler.dataset_validation import (
    HUMAN_REFERENCE_GENES,
    ECOLI_REFERENCE_GENES,
    YEAST_REFERENCE_GENES,
    SYNTHETIC_BENCHMARKS,
    format_dataset_report_text,
)
from biocompiler.benchmark import run_benchmarks, format_benchmark_report_text
from biocompiler.maxentscan import max_donor_score, max_acceptor_score, score_donor
from biocompiler.constants import AA_TO_CODONS, RESTRICTION_ENZYMES


def separator(title):
    print(f"\n{'='*80}")
    print(f"  {title}")
    print(f"{'='*80}\n")


def sub_separator(title):
    print(f"\n--- {title} ---\n")


# ============================================================================
# 1. FULL DATASET VALIDATION
# ============================================================================

separator("1. FULL DATASET VALIDATION (all organisms, all genes)")

print("Running comprehensive dataset validation...\n")
report = run_dataset_validation(
    include_cross_organism=True,
    include_optimization_improvement=True,
)
print(format_dataset_report_text(report))

# ============================================================================
# 2. BENCHMARK SUITE (HBB, INS, EGFP with real pre-mRNA)
# ============================================================================

separator("2. BENCHMARK SUITE (HBB, INS, EGFP)")

bench_report = run_benchmarks(include_optimization=True)
print(format_benchmark_report_text(bench_report))

# ============================================================================
# 3. DETAILED PER-GENE OPTIMIZATION ANALYSIS
# ============================================================================

separator("3. DETAILED PER-GENE OPTIMIZATION ANALYSIS")

all_genes = {}
for ds_name, ds in [("human", HUMAN_REFERENCE_GENES), ("ecoli", ECOLI_REFERENCE_GENES),
                     ("yeast", YEAST_REFERENCE_GENES), ("synthetic", SYNTHETIC_BENCHMARKS)]:
    for gene_name, gene_data in ds.items():
        all_genes[f"{ds_name}/{gene_name}"] = gene_data

print(f"Total genes to optimize: {len(all_genes)}\n")

results_table = []
for label, gene in sorted(all_genes.items()):
    protein = gene["protein"]
    organism = gene["organism"]
    t0 = time.perf_counter()

    # Optimize WITHOUT mutagenesis
    opt = optimize_sequence(
        target_protein=protein,
        organism=organism,
        gc_lo=0.30,
        gc_hi=0.70,
        cai_threshold=0.2,
        enable_mutagenesis=False,
    )
    elapsed = (time.perf_counter() - t0) * 1000

    # Count GT dinucleotides in the sequence
    gt_count = sum(1 for i in range(len(opt.sequence) - 1) if opt.sequence[i:i+2] == "GT")

    # Max splice scores
    max_d = max_donor_score(opt.sequence)
    max_a = max_acceptor_score(opt.sequence)

    # Count Valine positions
    valine_positions = [i for i, aa in enumerate(protein) if aa == "V"]

    results_table.append({
        "gene": label,
        "protein_len": len(protein),
        "organism": organism.replace("_", " "),
        "cai": opt.cai,
        "gc": opt.gc_content,
        "satisfied": len(opt.satisfied_predicates),
        "failed": len(opt.failed_predicates),
        "failed_list": opt.failed_predicates,
        "gt_count": gt_count,
        "valine_count": len(valine_positions),
        "max_donor": max_d,
        "max_acceptor": max_a,
        "fallback": opt.fallback_used,
        "time_ms": elapsed,
    })

# Print table
print(f"{'Gene':<30} {'Len':>4} {'Organism':<20} {'CAI':>6} {'GC':>5} {'Sat':>3} {'Fail':>4} "
      f"{'GTs':>3} {'Val':>3} {'MaxD':>6} {'MaxA':>6} {'ms':>7}")
print("-" * 120)

total_failures = 0
genes_with_failures = []

for r in sorted(results_table, key=lambda x: -x["failed"]):
    print(f"{r['gene']:<30} {r['protein_len']:>4} {r['organism']:<20} "
          f"{r['cai']:>6.4f} {r['gc']:>5.3f} {r['satisfied']:>3} {r['failed']:>4} "
          f"{r['gt_count']:>3} {r['valine_count']:>3} {r['max_donor']:>6.2f} "
          f"{r['max_acceptor']:>6.2f} {r['time_ms']:>7.1f}")
    if r["failed"] > 0:
        total_failures += r["failed"]
        genes_with_failures.append(r)
        print(f"  └─ Failed predicates: {r['failed_list']}")

print(f"\nGenes with predicate failures: {len(genes_with_failures)}/{len(all_genes)}")
print(f"Total failed predicates: {total_failures}")

# ============================================================================
# 4. VALENCE ANALYSIS — GT-MANDATORY AMINO ACIDS
# ============================================================================

separator("4. GT-MANDATORY AMINO ACID ANALYSIS")

print("Amino acids whose ALL codons contain the GT dinucleotide:")
for aa in sorted(GT_MANDATORY_AAS):
    codons = AA_TO_CODONS[aa]
    print(f"  {aa}: {', '.join(codons)} — ALL contain GT")
    for codon in codons:
        gt_pos = codon.find("GT")
        print(f"      {codon}: GT at position {gt_pos}")

print(f"\nAmino acids with GT-free codons available:")
for aa in sorted(set(AA_TO_CODONS.keys()) - GT_MANDATORY_AAS):
    codons = AA_TO_CODONS[aa]
    gt_free = [c for c in codons if "GT" not in c]
    gt_contain = [c for c in codons if "GT" in c]
    if gt_contain:
        print(f"  {aa}: GT-free={gt_free}, GT-containing={gt_contain}")

# ============================================================================
# 5. HBB DEEP DIVE — THE KEY MUTAGENESIS PROOF
# ============================================================================

separator("5. HBB DEEP DIVE — TYPE-DIRECTED MUTAGENESIS PROOF")

hbb_protein = HUMAN_REFERENCE_GENES["HBB"]["protein"]
hbb_organism = HUMAN_REFERENCE_GENES["HBB"]["organism"]

print(f"HBB protein ({len(hbb_protein)} aa):")
print(f"  {hbb_protein}")
print(f"  Valine positions: {[i for i, aa in enumerate(hbb_protein) if aa == 'V']}")
print(f"  Valine count: {sum(1 for aa in hbb_protein if aa == 'V')}/{len(hbb_protein)}")

# 5a. Optimize WITHOUT mutagenesis
sub_separator("5a. HBB optimization WITHOUT mutagenesis")
t0 = time.perf_counter()
hbb_opt_no_mut = optimize_sequence(
    target_protein=hbb_protein,
    organism=hbb_organism,
    gc_lo=0.30,
    gc_hi=0.70,
    cai_threshold=0.2,
    cryptic_splice_threshold=3.0,
    enable_mutagenesis=False,
)
elapsed_no_mut = (time.perf_counter() - t0) * 1000

print(f"  CAI: {hbb_opt_no_mut.cai:.4f}")
print(f"  GC:  {hbb_opt_no_mut.gc_content:.4f}")
print(f"  Satisfied predicates: {hbb_opt_no_mut.satisfied_predicates}")
print(f"  Failed predicates:    {hbb_opt_no_mut.failed_predicates}")
print(f"  Time: {elapsed_no_mut:.1f} ms")

# Find unrepairable cryptic donors
sub_separator("5b. Unrepairable cryptic splice donors in HBB")
unrepairable = find_unrepairable_cryptic_donors(
    hbb_opt_no_mut.sequence, hbb_protein, hbb_organism, threshold=3.0
)

print(f"  Unrepairable donor sites: {len(unrepairable)}")
for seq_pos, codon_idx, aa, score, fixable, gt_mandatory in unrepairable:
    codon = hbb_opt_no_mut.sequence[codon_idx*3:codon_idx*3+3]
    print(f"    Seq pos {seq_pos}, codon {codon_idx} ({aa}, codon={codon}), "
          f"MaxEntScan={score:.2f}, fixable_by_codon_swap={fixable}, gt_mandatory={gt_mandatory}")

# 5c. Optimize WITH mutagenesis
sub_separator("5c. HBB optimization WITH type-directed mutagenesis")
t0 = time.perf_counter()
hbb_opt_mut = optimize_sequence(
    target_protein=hbb_protein,
    organism=hbb_organism,
    gc_lo=0.30,
    gc_hi=0.70,
    cai_threshold=0.2,
    cryptic_splice_threshold=3.0,
    enable_mutagenesis=True,
    max_mutagenesis_substitutions=30,
    min_blosum62=-1,
)
elapsed_mut = (time.perf_counter() - t0) * 1000

print(f"  CAI: {hbb_opt_mut.cai:.4f}")
print(f"  GC:  {hbb_opt_mut.gc_content:.4f}")
print(f"  Satisfied predicates: {hbb_opt_mut.satisfied_predicates}")
print(f"  Failed predicates:    {hbb_opt_mut.failed_predicates}")
print(f"  Mutagenesis applied:  {hbb_opt_mut.mutagenesis_applied}")
if hbb_opt_mut.aa_substitutions:
    print(f"  Amino acid substitutions:")
    for sub in hbb_opt_mut.aa_substitutions:
        print(f"    Position {sub['position']}: {sub['from']}→{sub['to']} "
              f"(BLOSUM62={sub['blosum62']:+d}, reason: {sub['reason'][:80]})")
print(f"  Time: {elapsed_mut:.1f} ms")

# 5d. Comparison
sub_separator("5d. HBB Comparison: Without vs With Mutagenesis")
cai_cost = hbb_opt_no_mut.cai - hbb_opt_mut.cai
print(f"  {'Metric':<30} {'No Mutagenesis':>15} {'With Mutagenesis':>18}")
print(f"  {'-'*63}")
print(f"  {'CAI':<30} {hbb_opt_no_mut.cai:>15.4f} {hbb_opt_mut.cai:>18.4f}")
print(f"  {'GC':<30} {hbb_opt_no_mut.gc_content:>15.4f} {hbb_opt_mut.gc_content:>18.4f}")
print(f"  {'Failed predicates':<30} {len(hbb_opt_no_mut.failed_predicates):>15} {len(hbb_opt_mut.failed_predicates):>18}")
print(f"  {'Satisfied predicates':<30} {len(hbb_opt_no_mut.satisfied_predicates):>15} {len(hbb_opt_mut.satisfied_predicates):>18}")
if hbb_opt_mut.aa_substitutions:
    n_subs = len(hbb_opt_mut.aa_substitutions)
    identity = (1 - n_subs/len(hbb_protein)) * 100
    print(f"  {'AA substitutions':<30} {'0':>15} {n_subs:>18}")
    print(f"  {'Protein identity':<30} {'100.0%':>15} {identity:>17.1f}%")
    print(f"  {'CAI cost of mutagenesis':<30} {'N/A':>15} {cai_cost:>+18.4f}")

# ============================================================================
# 6. CROSS-ORGANISM ANALYSIS
# ============================================================================

separator("6. CROSS-ORGANISM OPTIMIZATION ANALYSIS")

test_proteins = {
    "HBB": HUMAN_REFERENCE_GENES["HBB"]["protein"],
    "GFP": ECOLI_REFERENCE_GENES["GFP"]["protein"],
}

organisms = ["Homo_sapiens", "Escherichia_coli", "Saccharomyces_cerevisiae"]

for gene_name, protein in test_proteins.items():
    sub_separator(f"{gene_name} optimized for different organisms")
    print(f"  {'Organism':<28} {'CAI':>6} {'GC':>5} {'Fail':>4} {'MaxDonor':>8} {'MaxAcceptor':>11}")
    print(f"  {'-'*62}")
    for org in organisms:
        opt = optimize_sequence(
            target_protein=protein,
            organism=org,
            gc_lo=0.20,
            gc_hi=0.80,
            cai_threshold=0.2,
            enable_mutagenesis=False,
        )
        md = max_donor_score(opt.sequence)
        ma = max_acceptor_score(opt.sequence)
        org_short = org.replace("Homo_sapiens", "Human").replace("Escherichia_coli", "E.coli").replace("Saccharomyces_cerevisiae", "Yeast")
        print(f"  {org_short:<28} {opt.cai:>6.4f} {opt.gc_content:>5.3f} "
              f"{len(opt.failed_predicates):>4} {md:>8.2f} {ma:>11.2f}")

# ============================================================================
# 7. RESTRICTION SITE ANALYSIS
# ============================================================================

separator("7. RESTRICTION SITE ELIMINATION ANALYSIS")

print("Enzymes in the restriction map:")
for name, site in sorted(RESTRICTION_ENZYMES.items()):
    print(f"  {name:<12} {site}")

print("\nOptimization success rate for restriction site elimination:")
for label, gene in sorted(all_genes.items()):
    protein = gene["protein"]
    organism = gene["organism"]
    opt = optimize_sequence(
        target_protein=protein,
        organism=organism,
        gc_lo=0.30,
        gc_hi=0.70,
        cai_threshold=0.2,
        enable_mutagenesis=False,
    )
    # Check if any restriction sites remain
    from biocompiler.constants import reverse_complement
    remaining = []
    for name, site in RESTRICTION_ENZYMES.items():
        site_upper = site.upper()
        rc = reverse_complement(site_upper)
        # Only check concrete sites
        if any(b not in "ACGT" for b in site_upper):
            continue
        if site_upper in opt.sequence or rc in opt.sequence:
            remaining.append(name)

    if remaining:
        print(f"  {label}: REMAINING sites: {remaining}")

# ============================================================================
# 8. MUTAGENESIS ACROSS ALL FAILING GENES
# ============================================================================

separator("8. MUTAGENESIS IMPACT ON ALL FAILING GENES")

for r in genes_with_failures:
    label = r["gene"]
    gene_data = all_genes[label]
    protein = gene_data["protein"]
    organism = gene_data["organism"]

    sub_separator(f"Mutagenesis for {label}")

    # Run with mutagenesis
    t0 = time.perf_counter()
    opt_mut = optimize_sequence(
        target_protein=protein,
        organism=organism,
        gc_lo=0.30,
        gc_hi=0.70,
        cai_threshold=0.2,
        cryptic_splice_threshold=3.0,
        enable_mutagenesis=True,
        max_mutagenesis_substitutions=30,
        min_blosum62=-1,
    )
    elapsed = (time.perf_counter() - t0) * 1000

    cai_before = r["cai"]
    cai_after = opt_mut.cai
    cai_delta = cai_after - cai_before

    print(f"  CAI:     {cai_before:.4f} → {cai_after:.4f} (Δ={cai_delta:+.4f})")
    print(f"  GC:      {r['gc']:.4f} → {opt_mut.gc_content:.4f}")
    print(f"  Failed:  {r['failed']} → {len(opt_mut.failed_predicates)}")
    print(f"  Applied: {opt_mut.mutagenesis_applied}")

    if opt_mut.aa_substitutions:
        # Group substitutions by type
        sub_types = {}
        for sub in opt_mut.aa_substitutions:
            key = f"{sub['from']}→{sub['to']}"
            sub_types[key] = sub_types.get(key, 0) + 1
        print(f"  Substitution types: {dict(sub_types)}")
        print(f"  Total substitutions: {len(opt_mut.aa_substitutions)}/{len(protein)} positions "
              f"({100*len(opt_mut.aa_substitutions)/len(protein):.1f}%)")
        print(f"  BLOSUM62 scores: {[s['blosum62'] for s in opt_mut.aa_substitutions]}")
        if opt_mut.aa_substitutions:
            avg_blosum = sum(s['blosum62'] for s in opt_mut.aa_substitutions) / len(opt_mut.aa_substitutions)
            print(f"  Average BLOSUM62: {avg_blosum:+.2f}")
    else:
        print(f"  No substitutions proposed (predicate may not be addressable by mutagenesis)")

    if opt_mut.failed_predicates:
        print(f"  Still failing: {opt_mut.failed_predicates}")

# ============================================================================
# 9. SUMMARY STATISTICS
# ============================================================================

separator("9. SUMMARY STATISTICS")

# Aggregate
all_cai = [r["cai"] for r in results_table]
all_gc = [r["gc"] for r in results_table]
all_fail = [r["failed"] for r in results_table]
all_time = [r["time_ms"] for r in results_table]

print(f"Genes tested:           {len(results_table)}")
print(f"Avg CAI:                {sum(all_cai)/len(all_cai):.4f}")
print(f"Min CAI:                {min(all_cai):.4f}")
print(f"Max CAI:                {max(all_cai):.4f}")
print(f"Avg GC:                 {sum(all_gc)/len(all_gc):.4f}")
print(f"Avg failed predicates:  {sum(all_fail)/len(all_fail):.2f}")
print(f"Genes with 0 failures:  {sum(1 for f in all_fail if f == 0)}/{len(all_fail)}")
print(f"Avg optimization time:  {sum(all_time)/len(all_time):.1f} ms")
print(f"Max optimization time:  {max(all_time):.1f} ms")

# By organism
sub_separator("By Organism")
for org_short in ["Human", "E.coli", "Yeast"]:
    org_results = [r for r in results_table if org_short.lower().replace(".", "") in r["organism"].lower()]
    if org_results:
        avg_cai = sum(r["cai"] for r in org_results) / len(org_results)
        avg_gc = sum(r["gc"] for r in org_results) / len(org_results)
        avg_fail = sum(r["failed"] for r in org_results) / len(org_results)
        print(f"  {org_short:<12} n={len(org_results):>2}, avg_CAI={avg_cai:.4f}, "
              f"avg_GC={avg_gc:.4f}, avg_failures={avg_fail:.2f}")

# Key insight
separator("KEY INSIGHT")
print("""
Valine is the ONLY amino acid whose ALL 4 codons (GTT, GTC, GTA, GTG)
contain the GT dinucleotide — the universal splice donor signal.

This means: when a Valine position creates an unrepairable cryptic splice
donor, NO codon swap can fix it. The ONLY solution is to change the
amino acid itself (type-directed mutagenesis).

V→I substitution (BLOSUM62 = +3) is the most conservative option:
- Both are branched-chain hydrophobic amino acids
- Similar volume (Val=140, Ile=166 Å³)
- Similar hydrophobicity (Kyte-Doolittle: Val=4.2, Ile=4.5)
- Isoleucine has GT-free codons: ATT, ATC, ATA

This is the central dogma crossing: DNA constraints (GT dinucleotide)
REQUIRE protein-level changes (V→I). The type predicate doesn't just
VERIFY — it DIRECTS the design.
""")

print("\n" + "="*80)
print("  ANALYSIS COMPLETE")
print("="*80)
