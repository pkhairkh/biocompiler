#!/usr/bin/env python3
"""
Full End-to-End Benchmark: BioCompiler vs DNAchisel vs Naive Baselines
=======================================================================
Comprehensive head-to-head comparison across:
  - 20 E. coli genes (prokaryotic optimization)
  - 5 human therapeutic genes (eukaryotic optimization)
  - Stress-test genes (extreme GC, long sequences)
  
Metrics: CAI (validated), GC%, runtime, constraint violations, codon diversity

CAI is computed using compute_cai_validated for BOTH tools — DNAchisel's
own CAI output is NOT trusted, ensuring fair head-to-head comparison.
"""

import sys
import time
import json
import statistics
import traceback
from pathlib import Path
from dataclasses import dataclass, field, asdict
from typing import Optional

# ── BioCompiler imports ──
from biocompiler.optimization import optimize_sequence, OptimizationResult
from biocompiler.translation import translate
from biocompiler.scanner import gc_content, validate_dna_sequence
from biocompiler.constants import AA_TO_CODONS, CODON_TABLE
from biocompiler.organisms import CODON_ADAPTIVENESS_TABLES, PREFERRED_CODON_TABLES, SUPPORTED_ORGANISMS
from biocompiler.benchmarking.metrics import compute_cai_validated

# ── DNAchisel imports ──
import dnachisel as dc

# ── Gene sets ──
from biocompiler.benchmarking.gene_sets import (
    E_COLI_EXTENDED, HUMAN_THERAPEUTIC, VACCINE_ANTIGEN_GENES, STRESS_TEST_GENES
)


@dataclass
class BenchmarkResult:
    tool: str
    gene: str
    organism: str
    cais: dict = field(default_factory=dict)   # organism -> CAI
    gc_content: float = 0.0
    runtime_s: float = 0.0
    dna_length: int = 0
    constraint_violations: list = field(default_factory=list)
    error: str = ""
    sequence: str = ""


def _organism_key(name: str) -> str:
    """Normalize organism name to BioCompiler key."""
    mapping = {
        "Escherichia coli": "Escherichia_coli",
        "Homo sapiens": "Homo_sapiens",
        "Saccharomyces cerevisiae": "Saccharomyces_cerevisiae",
        "Mus musculus": "Mus_musculus",
        "CHO K1": "CHO_K1",
    }
    return mapping.get(name, name.replace(" ", "_"))


def _dnachisel_species(organism_key: str) -> str:
    """Map BioCompiler organism key to DNAchisel species name."""
    mapping = {
        "Escherichia_coli": "e_coli",
        "Homo_sapiens": "h_sapiens",
        "Saccharomyces_cerevisiae": "s_cerevisiae",
        "Mus_musculus": "m_musculus",
    }
    return mapping.get(organism_key, organism_key.lower())


# ═══════════════════════════════════════════════════════════════
# Tool 1: Naive baseline — just pick most-frequent codon per AA
# ═══════════════════════════════════════════════════════════════

def optimize_naive(protein: str, organism: str) -> str:
    """Return DNA using most-frequent codon for each amino acid."""
    preferred = PREFERRED_CODON_TABLES.get(organism, {})
    if not preferred:
        # fallback: first codon per AA
        return "".join(AA_TO_CODONS[aa][0] for aa in protein)
    seq_parts = []
    for aa in protein:
        codon = preferred.get(aa, AA_TO_CODONS[aa][0])
        seq_parts.append(codon)
    return "".join(seq_parts)


# ═══════════════════════════════════════════════════════════════
# Tool 2: BioCompiler (greedy optimizer)
# ═══════════════════════════════════════════════════════════════

def optimize_biocompiler(protein: str, organism: str, gc_lo: float = 0.30, gc_hi: float = 0.70) -> BenchmarkResult:
    result = BenchmarkResult(tool="biocompiler", gene="", organism=organism)
    t0 = time.perf_counter()
    try:
        opt_result = optimize_sequence(
            target_protein=protein,
            organism=organism,
            gc_lo=gc_lo,
            gc_hi=gc_hi,
        )
        t1 = time.perf_counter()
        result.runtime_s = t1 - t0
        result.sequence = opt_result.sequence
        result.dna_length = len(opt_result.sequence)
        result.gc_content = gc_content(opt_result.sequence)
        result.constraint_violations = list(opt_result.failed_predicates) if opt_result.failed_predicates else []
        # Compute CAI for all supported organisms using validated evaluator
        for org in SUPPORTED_ORGANISMS:
            try:
                result.cais[org] = compute_cai_validated(opt_result.sequence, org)
            except Exception:
                pass
    except Exception as e:
        t1 = time.perf_counter()
        result.runtime_s = t1 - t0
        result.error = f"{type(e).__name__}: {e}"
    return result


# ═══════════════════════════════════════════════════════════════
# Tool 3: DNAchisel
# ═══════════════════════════════════════════════════════════════

def optimize_dnachisel(protein: str, organism: str, gc_lo: float = 0.30, gc_hi: float = 0.70) -> BenchmarkResult:
    result = BenchmarkResult(tool="dnachisel", gene="", organism=organism)
    species = _dnachisel_species(organism)
    
    # Build starting sequence
    preferred = PREFERRED_CODON_TABLES.get(organism, {})
    if preferred:
        start_seq = "".join(preferred.get(aa, AA_TO_CODONS[aa][0]) for aa in protein)
    else:
        start_seq = "".join(AA_TO_CODONS[aa][0] for aa in protein)
    
    t0 = time.perf_counter()
    try:
        constraints = [
            dc.EnforceTranslation(translation=protein),
            dc.EnforceGCContent(mini=gc_lo, maxi=gc_hi),
        ]
        objectives = [
            dc.CodonOptimize(species=species),
        ]
        
        problem = dc.DnaOptimizationProblem(
            sequence=start_seq,
            constraints=constraints,
            objectives=objectives,
            logger=None,
        )
        problem.resolve_constraints()
        problem.optimize()
        t1 = time.perf_counter()
        
        result.runtime_s = t1 - t0
        result.sequence = problem.sequence
        result.dna_length = len(problem.sequence)
        result.gc_content = gc_content(problem.sequence)
        
        # Compute CAI using validated evaluator (fair comparison)
        # DNAchisel's own CAI output is NOT trusted
        for org in SUPPORTED_ORGANISMS:
            try:
                result.cais[org] = compute_cai_validated(problem.sequence, org)
            except Exception:
                pass
    except Exception as e:
        t1 = time.perf_counter()
        result.runtime_s = t1 - t0
        result.error = f"{type(e).__name__}: {e}"
    return result


# ═══════════════════════════════════════════════════════════════
# Tool 4: DNAchisel with restriction sites
# ═══════════════════════════════════════════════════════════════

def optimize_dnachisel_with_re(protein: str, organism: str, gc_lo: float = 0.30, gc_hi: float = 0.70) -> BenchmarkResult:
    """DNAchisel with common restriction enzyme avoidance (EcoRI, BamHI, XhoI, NotI)."""
    result = BenchmarkResult(tool="dnachisel+re", gene="", organism=organism)
    species = _dnachisel_species(organism)
    
    preferred = PREFERRED_CODON_TABLES.get(organism, {})
    if preferred:
        start_seq = "".join(preferred.get(aa, AA_TO_CODONS[aa][0]) for aa in protein)
    else:
        start_seq = "".join(AA_TO_CODONS[aa][0] for aa in protein)
    
    t0 = time.perf_counter()
    try:
        constraints = [
            dc.EnforceTranslation(translation=protein),
            dc.EnforceGCContent(mini=gc_lo, maxi=gc_hi),
            dc.AvoidPattern("GAATTC"),  # EcoRI
            dc.AvoidPattern("GGATCC"),  # BamHI
            dc.AvoidPattern("CTCGAG"),  # XhoI
            dc.AvoidPattern("GCGGCCGC"),  # NotI
        ]
        objectives = [
            dc.CodonOptimize(species=species),
        ]
        
        problem = dc.DnaOptimizationProblem(
            sequence=start_seq,
            constraints=constraints,
            objectives=objectives,
            logger=None,
        )
        problem.resolve_constraints()
        problem.optimize()
        t1 = time.perf_counter()
        
        result.runtime_s = t1 - t0
        result.sequence = problem.sequence
        result.dna_length = len(problem.sequence)
        result.gc_content = gc_content(problem.sequence)
        
        # Compute CAI using validated evaluator (fair comparison)
        # DNAchisel's own CAI output is NOT trusted
        for org in SUPPORTED_ORGANISMS:
            try:
                result.cais[org] = compute_cai_validated(problem.sequence, org)
            except Exception:
                pass
    except Exception as e:
        t1 = time.perf_counter()
        result.runtime_s = t1 - t0
        result.error = f"{type(e).__name__}: {e}"
    return result


# ═══════════════════════════════════════════════════════════════
# Tool 5: BioCompiler with organism-aware (disable splice/CpG for prokaryotes)
# ═══════════════════════════════════════════════════════════════

def optimize_biocompiler_organism_aware(protein: str, organism: str, gc_lo: float = 0.30, gc_hi: float = 0.70) -> BenchmarkResult:
    """BioCompiler with organism-aware constraint selection."""
    result = BenchmarkResult(tool="biocompiler-aware", gene="", organism=organism)
    t0 = time.perf_counter()
    try:
        # Disable eukaryotic-only constraints for prokaryotes
        is_eukaryote = organism in ("Homo_sapiens", "Mus_musculus", "CHO_K1", "Saccharomyces_cerevisiae")
        
        # Use lower-level API to control which constraints apply
        from biocompiler.optimization import BioOptimizer
        opt = BioOptimizer(
            species=organism.replace("_", " ").lower() if is_eukaryote else organism.replace("_", " ").lower(),
            enzymes=[],  # skip restriction sites for fair comparison
            strategy="cai_first" if not is_eukaryote else "constraint_first",
        )
        # For prokaryotes: skip CpG and splice constraints
        dna_start = "".join(PREFERRED_CODON_TABLES.get(organism, {}).get(aa, AA_TO_CODONS[aa][0]) for aa in protein) + "TAA"
        optimized, pred_results, cert = opt.optimize(dna_start)
        t1 = time.perf_counter()
        
        result.runtime_s = t1 - t0
        result.sequence = optimized
        result.dna_length = len(optimized)
        result.gc_content = gc_content(optimized)
        
        for org in SUPPORTED_ORGANISMS:
            try:
                result.cais[org] = compute_cai_validated(optimized, org)
            except Exception:
                pass
    except Exception as e:
        t1 = time.perf_counter()
        result.runtime_s = t1 - t0
        result.error = f"{type(e).__name__}: {e}"
    return result


# ═══════════════════════════════════════════════════════════════
# Main benchmark runner
# ═══════════════════════════════════════════════════════════════

def run_full_benchmark():
    all_results: list[BenchmarkResult] = []
    
    # ── 1. E. coli genes (prokaryotic) ──
    print("\n" + "="*80)
    print("SECTION 1: E. coli Prokaryotic Optimization (20 genes)")
    print("="*80)
    
    for gene_name, gene_data in E_COLI_EXTENDED.items():
        protein = gene_data["protein"]
        org_key = _organism_key(gene_data.get("organism", "Escherichia coli"))
        print(f"\n  Gene: {gene_name} ({len(protein)} aa)")
        
        # Naive
        r = BenchmarkResult(tool="naive", gene=gene_name, organism=org_key)
        t0 = time.perf_counter()
        r.sequence = optimize_naive(protein, org_key)
        r.runtime_s = time.perf_counter() - t0
        r.dna_length = len(r.sequence)
        r.gc_content = gc_content(r.sequence)
        for o in SUPPORTED_ORGANISMS:
            try: r.cais[o] = compute_cai_validated(r.sequence, o)
            except: pass
        all_results.append(r)
        print(f"    naive:          CAI={r.cais.get(org_key, 0):.4f}  GC={r.gc_content:.4f}  t={r.runtime_s*1000:.1f}ms")
        
        # BioCompiler
        r = optimize_biocompiler(protein, org_key)
        r.gene = gene_name
        all_results.append(r)
        print(f"    biocompiler:    CAI={r.cais.get(org_key, 0):.4f}  GC={r.gc_content:.4f}  t={r.runtime_s*1000:.1f}ms  violations={r.constraint_violations}{f'  ERROR: {r.error[:60]}' if r.error else ''}")
        
        # DNAchisel
        r = optimize_dnachisel(protein, org_key)
        r.gene = gene_name
        all_results.append(r)
        print(f"    dnachisel:      CAI={r.cais.get(org_key, 0):.4f}  GC={r.gc_content:.4f}  t={r.runtime_s*1000:.1f}ms{f'  ERROR: {r.error[:60]}' if r.error else ''}")
        
        # DNAchisel + RE
        r = optimize_dnachisel_with_re(protein, org_key)
        r.gene = gene_name
        all_results.append(r)
        print(f"    dnachisel+re:   CAI={r.cais.get(org_key, 0):.4f}  GC={r.gc_content:.4f}  t={r.runtime_s*1000:.1f}ms{f'  ERROR: {r.error[:60]}' if r.error else ''}")
    
    # ── 2. Human therapeutic genes (eukaryotic) ──
    print("\n" + "="*80)
    print("SECTION 2: Human Therapeutic Eukaryotic Optimization (5 genes)")
    print("="*80)
    
    for gene_name, gene_data in HUMAN_THERAPEUTIC.items():
        protein = gene_data["protein"]
        org_key = _organism_key(gene_data.get("organism", "Homo sapiens"))
        print(f"\n  Gene: {gene_name} ({len(protein)} aa)")
        
        # Naive
        r = BenchmarkResult(tool="naive", gene=gene_name, organism=org_key)
        t0 = time.perf_counter()
        r.sequence = optimize_naive(protein, org_key)
        r.runtime_s = time.perf_counter() - t0
        r.dna_length = len(r.sequence)
        r.gc_content = gc_content(r.sequence)
        for o in SUPPORTED_ORGANISMS:
            try: r.cais[o] = compute_cai_validated(r.sequence, o)
            except: pass
        all_results.append(r)
        print(f"    naive:          CAI={r.cais.get(org_key, 0):.4f}  GC={r.gc_content:.4f}")
        
        # BioCompiler
        r = optimize_biocompiler(protein, org_key)
        r.gene = gene_name
        all_results.append(r)
        print(f"    biocompiler:    CAI={r.cais.get(org_key, 0):.4f}  GC={r.gc_content:.4f}  t={r.runtime_s*1000:.1f}ms  violations={r.constraint_violations}{f'  ERROR: {r.error[:60]}' if r.error else ''}")
        
        # DNAchisel
        r = optimize_dnachisel(protein, org_key)
        r.gene = gene_name
        all_results.append(r)
        print(f"    dnachisel:      CAI={r.cais.get(org_key, 0):.4f}  GC={r.gc_content:.4f}  t={r.runtime_s*1000:.1f}ms{f'  ERROR: {r.error[:60]}' if r.error else ''}")
    
    # ── 3. Vaccine antigens (designed for human expression) ──
    print("\n" + "="*80)
    print("SECTION 3: Vaccine Antigen Optimization (selected genes, human expression)")
    print("="*80)
    
    vaccine_subset = {"SARS2_RBD", "RABV_G", "Zika_E"}
    for gene_name, gene_data in VACCINE_ANTIGEN_GENES.items():
        if gene_name not in vaccine_subset:
            continue
        protein = gene_data["protein_sequence"]
        org_key = "Homo_sapiens"  # target human expression
        print(f"\n  Gene: {gene_name} ({len(protein)} aa)")
        
        # Naive
        r = BenchmarkResult(tool="naive", gene=gene_name, organism=org_key)
        t0 = time.perf_counter()
        r.sequence = optimize_naive(protein, org_key)
        r.runtime_s = time.perf_counter() - t0
        r.dna_length = len(r.sequence)
        r.gc_content = gc_content(r.sequence)
        for o in SUPPORTED_ORGANISMS:
            try: r.cais[o] = compute_cai_validated(r.sequence, o)
            except: pass
        all_results.append(r)
        print(f"    naive:          CAI={r.cais.get(org_key, 0):.4f}  GC={r.gc_content:.4f}")
        
        # BioCompiler
        r = optimize_biocompiler(protein, org_key)
        r.gene = gene_name
        all_results.append(r)
        print(f"    biocompiler:    CAI={r.cais.get(org_key, 0):.4f}  GC={r.gc_content:.4f}  t={r.runtime_s*1000:.1f}ms  violations={r.constraint_violations}{f'  ERROR: {r.error[:60]}' if r.error else ''}")
        
        # DNAchisel
        r = optimize_dnachisel(protein, org_key)
        r.gene = gene_name
        all_results.append(r)
        print(f"    dnachisel:      CAI={r.cais.get(org_key, 0):.4f}  GC={r.gc_content:.4f}  t={r.runtime_s*1000:.1f}ms{f'  ERROR: {r.error[:60]}' if r.error else ''}")
    
    # ── 4. Cross-organism: optimize EGFP for all organisms ──
    print("\n" + "="*80)
    print("SECTION 4: Cross-Organism EGFP Optimization")
    print("="*80)
    
    egfp = "MVSKGEELFTGVVPILVELDGDVNGHKFSVSGEGEGDATYGKLTLKFICTTGKLPVPWPTLVTTFSYGVQCFSRYPDHMKQHDFFKSAMPEGYVQERTIFFKDDGNYKTRAEVKFEGDTLVNRIELKGIDFKEDGNILGHKLEYNYNSHNVYITADKQKNGIKANFKIRHNIEDGSVQLADHYQQNTPIGDGPVLLPDNHYLSTQSALSKDPNEKRDHMVLLEFVTAAGITHGMDELYK"
    
    for org in SUPPORTED_ORGANISMS:
        print(f"\n  Organism: {org}")
        
        # Naive
        r = BenchmarkResult(tool="naive", gene="EGFP", organism=org)
        t0 = time.perf_counter()
        r.sequence = optimize_naive(egfp, org)
        r.runtime_s = time.perf_counter() - t0
        r.dna_length = len(r.sequence)
        r.gc_content = gc_content(r.sequence)
        for o in SUPPORTED_ORGANISMS:
            try: r.cais[o] = compute_cai_validated(r.sequence, o)
            except: pass
        all_results.append(r)
        print(f"    naive:          CAI={r.cais.get(org, 0):.4f}  GC={r.gc_content:.4f}")
        
        # BioCompiler
        r = optimize_biocompiler(egfp, org)
        r.gene = "EGFP"
        all_results.append(r)
        print(f"    biocompiler:    CAI={r.cais.get(org, 0):.4f}  GC={r.gc_content:.4f}  t={r.runtime_s*1000:.1f}ms  violations={r.constraint_violations}{f'  ERROR: {r.error[:60]}' if r.error else ''}")
        
        # DNAchisel
        r = optimize_dnachisel(egfp, org)
        r.gene = "EGFP"
        all_results.append(r)
        print(f"    dnachisel:      CAI={r.cais.get(org, 0):.4f}  GC={r.gc_content:.4f}  t={r.runtime_s*1000:.1f}ms{f'  ERROR: {r.error[:60]}' if r.error else ''}")
    
    # ── 5. Stress tests ──
    print("\n" + "="*80)
    print("SECTION 5: Stress-Test Genes")
    print("="*80)
    
    for gene_name, gene_data in STRESS_TEST_GENES.items():
        protein = gene_data["protein_sequence"]
        org_key = "Escherichia_coli"  # default to E. coli for stress tests
        print(f"\n  Gene: {gene_name} ({len(protein)} aa, category: {gene_data.get('stress_category', 'unknown')})")
        
        # BioCompiler
        r = optimize_biocompiler(protein, org_key)
        r.gene = gene_name
        all_results.append(r)
        print(f"    biocompiler:    CAI={r.cais.get(org_key, 0):.4f}  GC={r.gc_content:.4f}  t={r.runtime_s*1000:.1f}ms{f'  ERROR: {r.error[:60]}' if r.error else ''}")
        
        # DNAchisel
        r = optimize_dnachisel(protein, org_key)
        r.gene = gene_name
        all_results.append(r)
        print(f"    dnachisel:      CAI={r.cais.get(org_key, 0):.4f}  GC={r.gc_content:.4f}  t={r.runtime_s*1000:.1f}ms{f'  ERROR: {r.error[:60]}' if r.error else ''}")
    
    # ═══════════════════════════════════════════════════════════
    # Summary
    # ═══════════════════════════════════════════════════════════
    print("\n" + "="*80)
    print("SUMMARY")
    print("="*80)
    
    # Group results by tool and compute statistics
    from collections import defaultdict
    tool_results = defaultdict(list)
    for r in all_results:
        if not r.error:
            tool_results[r.tool].append(r)
    
    for tool, results in sorted(tool_results.items()):
        cais = [r.cais.get(r.organism, 0) for r in results]
        gcs = [r.gc_content for r in results]
        times = [r.runtime_s for r in results]
        
        print(f"\n  {tool}:")
        print(f"    Genes tested: {len(results)}")
        if cais:
            print(f"    CAI: mean={statistics.mean(cais):.4f}  median={statistics.median(cais):.4f}  min={min(cais):.4f}  max={max(cais):.4f}  stdev={statistics.stdev(cais):.4f}" if len(cais) > 1 else f"    CAI: {cais[0]:.4f}")
        if gcs:
            print(f"    GC%: mean={statistics.mean(gcs):.4f}  median={statistics.median(gcs):.4f}")
        if times:
            print(f"    Time: mean={statistics.mean(times)*1000:.1f}ms  median={statistics.median(times)*1000:.1f}ms  total={sum(times):.2f}s")
    
    # ── Head-to-head comparison per gene ──
    print("\n" + "-"*80)
    print("HEAD-TO-HEAD: BioCompiler vs DNAchisel (CAI on target organism)")
    print("-"*80)
    
    gene_tools = defaultdict(dict)
    for r in all_results:
        if not r.error and r.tool in ("biocompiler", "dnachisel"):
            gene_tools[r.gene][r.tool] = r.cais.get(r.organism, 0)
    
    bc_wins = 0
    dc_wins = 0
    ties = 0
    cai_deltas = []
    
    for gene, tools in sorted(gene_tools.items()):
        bc_cai = tools.get("biocompiler", 0)
        dc_cai = tools.get("dnachisel", 0)
        delta = bc_cai - dc_cai
        cai_deltas.append(delta)
        winner = "BC" if delta > 0.001 else ("DC" if delta < -0.001 else "TIE")
        if winner == "BC": bc_wins += 1
        elif winner == "DC": dc_wins += 1
        else: ties += 1
        print(f"    {gene:20s}  BC={bc_cai:.4f}  DC={dc_cai:.4f}  delta={delta:+.4f}  [{winner}]")
    
    print(f"\n  BioCompiler wins: {bc_wins}  DNAchisel wins: {dc_wins}  Ties: {ties}")
    if cai_deltas:
        print(f"  Mean CAI delta: {statistics.mean(cai_deltas):+.4f}")
        print(f"  Median CAI delta: {statistics.median(cai_deltas):+.4f}")
    
    # ── E. coli specific comparison ──
    print("\n" + "-"*80)
    print("E. COLI SPECIFIC: BioCompiler vs DNAchisel vs Naive")
    print("-"*80)
    
    ecoli_results = [r for r in all_results if r.organism == "Escherichia_coli" and not r.error]
    for tool in ["naive", "biocompiler", "dnachisel", "dnachisel+re"]:
        tool_ecoli = [r for r in ecoli_results if r.tool == tool]
        if tool_ecoli:
            cais = [r.cais.get("Escherichia_coli", 0) for r in tool_ecoli]
            print(f"  {tool:18s}: CAI mean={statistics.mean(cais):.4f}  range=[{min(cais):.4f}, {max(cais):.4f}]  n={len(cais)}")
    
    # ── Human specific comparison ──
    print("\n" + "-"*80)
    print("HUMAN SPECIFIC: BioCompiler vs DNAchisel vs Naive")
    print("-"*80)
    
    human_results = [r for r in all_results if r.organism == "Homo_sapiens" and not r.error]
    for tool in ["naive", "biocompiler", "dnachisel"]:
        tool_human = [r for r in human_results if r.tool == tool]
        if tool_human:
            cais = [r.cais.get("Homo_sapiens", 0) for r in tool_human]
            print(f"  {tool:18s}: CAI mean={statistics.mean(cais):.4f}  range=[{min(cais):.4f}, {max(cais):.4f}]  n={len(cais)}")
    
    # ── Speed comparison ──
    print("\n" + "-"*80)
    print("SPEED COMPARISON (all genes)")
    print("-"*80)
    
    for tool in ["naive", "biocompiler", "dnachisel", "dnachisel+re"]:
        tool_all = [r for r in all_results if r.tool == tool and not r.error]
        if tool_all:
            times = [r.runtime_s for r in tool_all]
            print(f"  {tool:18s}: total={sum(times):.2f}s  mean={statistics.mean(times)*1000:.1f}ms  median={statistics.median(times)*1000:.1f}ms")
    
    # ── Constraint violations ──
    print("\n" + "-"*80)
    print("CONSTRAINT VIOLATIONS (BioCompiler)")
    print("-"*80)
    
    bc_results = [r for r in all_results if r.tool == "biocompiler" and not r.error]
    violation_counts = defaultdict(int)
    for r in bc_results:
        for v in r.constraint_violations:
            violation_counts[v] += 1
    for v, count in sorted(violation_counts.items(), key=lambda x: -x[1]):
        print(f"  {v}: {count}/{len(bc_results)} genes")
    
    # ── Save JSON results ──
    output_dir = Path("/home/z/my-project/biocompiler/benchmark_results")
    output_dir.mkdir(exist_ok=True)
    
    json_data = []
    for r in all_results:
        d = asdict(r)
        d["cais"] = {k: round(v, 6) for k, v in d["cais"].items()}
        d["gc_content"] = round(d["gc_content"], 6)
        d["runtime_s"] = round(d["runtime_s"], 6)
        json_data.append(d)
    
    with open(output_dir / "e2e_benchmark_results.json", "w") as f:
        json.dump(json_data, f, indent=2)
    print(f"\n  Results saved to {output_dir / 'e2e_benchmark_results.json'}")
    
    return all_results


if __name__ == "__main__":
    results = run_full_benchmark()
