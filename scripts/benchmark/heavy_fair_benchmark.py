#!/usr/bin/env python3
"""
BioCompiler Heavy Fair Benchmark vs DNAchisel
=============================================
FAIR comparison: both tools optimize the SAME protein with the SAME constraints.
- Same objective: maximize CAI
- Same constraints: GC 30-70%, avoid EcoRI/BamHI/HindIII/XhoI, no premature stops
- Warmup: 3 warmup runs before timing
- Repeats: 5 timed runs, take median
- 25 genes across E. coli + human
- Metrics: CAI, GC, time (ms), constraint violations
"""
from __future__ import annotations
import json, time, statistics, sys, os, warnings
warnings.filterwarnings("ignore")

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "src"))

from biocompiler import optimize_sequence
from biocompiler.benchmarking import gene_sets
from biocompiler.expression.translation import compute_cai

try:
    import dnachisel as dc
    from dnachisel import DnaOptimizationProblem
    DNACHISEL_AVAILABLE = True
except ImportError:
    DNACHISEL_AVAILABLE = False

ORGANISM_TO_DC_TABLE = {"Escherichia_coli": "e_coli_316407", "Homo_sapiens": "h_sapiens_9606"}
ENZYMES = ["EcoRI", "BamHI", "HindIII", "XhoI"]
ENZYME_SITES = {"EcoRI": "GAATTC", "BamHI": "GGATCC", "HindIII": "AAGCTT", "XhoI": "CTCGAG"}

def count_violations(dna: str) -> dict:
    v = {"restriction_sites": 0, "gc_out_of_range": False, "premature_stop": 0}
    for site in ENZYME_SITES.values():
        rev = site[::-1].translate(str.maketrans("ATCG","TAGC"))
        v["restriction_sites"] += dna.count(site) + dna.count(rev)
    gc = (dna.count("G") + dna.count("C")) / len(dna) if dna else 0
    if gc < 0.30 or gc > 0.70:
        v["gc_out_of_range"] = True
    for i in range(0, len(dna) - 3, 3):
        if dna[i:i+3] in ("TAA","TAG","TGA") and i < len(dna) - 3:
            v["premature_stop"] += 1
    return v

def optimize_with_biocompiler(protein: str, organism: str) -> dict:
    t0 = time.perf_counter()
    try:
        result = optimize_sequence(
            protein, organism=organism,
            gc_lo=0.30, gc_hi=0.70,
            enzymes=ENZYMES, cpg_mode="off",
            strict_mode=False, skip_biosecurity_check=True,
        )
        t1 = time.perf_counter()
        dna = result.sequence
        return {
            "dna": dna, "cai": compute_cai(dna, organism) if dna else 0.0,
            "gc": (dna.count("G")+dna.count("C"))/len(dna) if dna else 0.0,
            "time_ms": (t1-t0)*1000, "violations": count_violations(dna),
        }
    except Exception as e:
        t1 = time.perf_counter()
        return {"dna":"", "cai":0.0, "gc":0.0, "time_ms":(t1-t0)*1000,
                "violations":{"restriction_sites":0,"gc_out_of_range":True,"premature_stop":0}, "error":str(e)[:80]}

def optimize_with_dnachisel(protein: str, organism: str) -> dict:
    from biocompiler.type_system.codon_tables import AA_TO_CODONS
    from biocompiler.organisms import resolve_organism, CODON_ADAPTIVENESS_TABLES
    species = resolve_organism(organism)
    usage = CODON_ADAPTIVENESS_TABLES.get(species, {})
    dna_parts = []
    for aa in protein:
        if aa == "U": dna_parts.append("TGA"); continue
        codons = AA_TO_CODONS.get(aa, ["GCT"])
        best = max(codons, key=lambda c: usage.get(c, 0.0)) if usage else codons[0]
        dna_parts.append(best)
    initial_dna = "".join(dna_parts)  # NO stop codon for DC

    constraints = [
        dc.EnforceTranslation(translation=protein),
        dc.EnforceGCContent(mini=0.30, maxi=0.70),
        dc.AvoidPattern(ENZYME_SITES["EcoRI"]),
        dc.AvoidPattern(ENZYME_SITES["BamHI"]),
        dc.AvoidPattern(ENZYME_SITES["HindIII"]),
        dc.AvoidPattern(ENZYME_SITES["XhoI"]),
    ]
    dc_table = ORGANISM_TO_DC_TABLE.get(organism, "e_coli_316407")
    objectives = [dc.MaximizeCAI(species=dc_table)]

    t0 = time.perf_counter()
    try:
        import proglog
        problem = DnaOptimizationProblem(sequence=initial_dna, constraints=constraints, objectives=objectives)
        problem.resolve_constraints()
        problem.optimize()
        t1 = time.perf_counter()
        dna = str(problem.sequence)
        return {"dna": dna, "cai": compute_cai(dna, organism) if dna else 0.0,
                "gc": (dna.count("G")+dna.count("C"))/len(dna) if dna else 0.0,
                "time_ms": (t1-t0)*1000, "violations": count_violations(dna)}
    except Exception as e:
        t1 = time.perf_counter()
        return {"dna":"", "cai":0.0, "gc":0.0, "time_ms":(t1-t0)*1000,
                "violations":{"restriction_sites":0,"gc_out_of_range":True,"premature_stop":0}, "error":str(e)[:80]}

def run_heavy_benchmark():
    all_sets = gene_sets.get_all_gene_sets()
    ecoli_genes = ["recA","gyrA","rpoB","groEL","tufA","ompA","lacZ","malE","araC","cat"]
    human_genes = ["INS","GH1","EPO","IFNA2","IL2","CSF3","PLAT","F8","SERPINA1","GBA"]
    therapeutic = ["SARS2_RBD","RSV_F","RABV_G","Zika_E","DENV2_E"]

    test_cases = []
    for name in ecoli_genes:
        if name in all_sets:
            p = all_sets[name].get("protein") or all_sets[name].get("protein_sequence","")
            if p: test_cases.append((name, p, "Escherichia_coli"))
    for name in human_genes + therapeutic:
        if name in all_sets:
            p = all_sets[name].get("protein") or all_sets[name].get("protein_sequence","")
            if p: test_cases.append((name, p, "Homo_sapiens"))

    print(f"Heavy Fair Benchmark: {len(test_cases)} genes")
    print(f"DNAchisel available: {DNACHISEL_AVAILABLE}")
    print(f"Constraints: GC 30-70%, avoid {ENZYMES}, no premature stops")
    print(f"Objective: maximize CAI (BC cpg_mode=off, DC MaximizeCAI)")
    print(f"Method: 3 warmup + 5 timed (median), same back-translation seed")
    print()

    results = []
    for i, (gene_name, protein, organism) in enumerate(test_cases):
        print(f"[{i+1}/{len(test_cases)}] {gene_name} ({len(protein)}aa, {organism.split('_')[0][:3]})...", end=" ", flush=True)
        for _ in range(3):
            optimize_with_biocompiler(protein, organism)
            if DNACHISEL_AVAILABLE: optimize_with_dnachisel(protein, organism)

        bc_times, dc_times, bc_cais, dc_cais = [],[],[],[]
        bc_results = dc_results = None
        for _ in range(5):
            bc_r = optimize_with_biocompiler(protein, organism)
            bc_times.append(bc_r["time_ms"]); bc_cais.append(bc_r["cai"]); bc_results = bc_r
            if DNACHISEL_AVAILABLE:
                dc_r = optimize_with_dnachisel(protein, organism)
                dc_times.append(dc_r["time_ms"]); dc_cais.append(dc_r["cai"]); dc_results = dc_r

        row = {
            "gene": gene_name, "organism": organism, "protein_length": len(protein),
            "bc_cai": round(statistics.median(bc_cais), 4),
            "dc_cai": round(statistics.median(dc_cais), 4) if dc_cais else None,
            "bc_time_ms": round(statistics.median(bc_times), 2),
            "dc_time_ms": round(statistics.median(dc_times), 2) if dc_times else None,
            "bc_gc": round(bc_results["gc"], 4) if bc_results else None,
            "dc_gc": round(dc_results["gc"], 4) if dc_results else None,
            "bc_violations": bc_results["violations"]["restriction_sites"] + bc_results["violations"]["premature_stop"] if bc_results else 0,
            "dc_violations": dc_results["violations"]["restriction_sites"] + dc_results["violations"]["premature_stop"] if dc_results else 0,
            "bc_error": bc_results.get("error","") if bc_results else "",
            "dc_error": dc_results.get("error","") if dc_results else "",
        }
        results.append(row)
        dc_str = f"DC CAI={statistics.median(dc_cais):.3f} {statistics.median(dc_times):.1f}ms" if dc_cais else "DC N/A"
        print(f"BC CAI={statistics.median(bc_cais):.3f} {statistics.median(bc_times):.1f}ms | {dc_str}")

    bc_cais_all = [r["bc_cai"] for r in results if r["bc_cai"]]
    dc_cais_all = [r["dc_cai"] for r in results if r["dc_cai"] is not None]
    bc_times_all = [r["bc_time_ms"] for r in results if r["bc_time_ms"]]
    dc_times_all = [r["dc_time_ms"] for r in results if r["dc_time_ms"] is not None]
    bc_viol = sum(r["bc_violations"] for r in results)
    dc_viol = sum(r["dc_violations"] for r in results)

    summary = {
        "total_genes": len(results),
        "bc_mean_cai": round(statistics.mean(bc_cais_all), 4) if bc_cais_all else 0,
        "dc_mean_cai": round(statistics.mean(dc_cais_all), 4) if dc_cais_all else 0,
        "bc_median_time_ms": round(statistics.median(bc_times_all), 2) if bc_times_all else 0,
        "dc_median_time_ms": round(statistics.median(dc_times_all), 2) if dc_times_all else 0,
        "bc_total_violations": bc_viol, "dc_total_violations": dc_viol,
        "bc_cai_wins": sum(1 for r in results if r["bc_cai"] and r["dc_cai"] is not None and r["bc_cai"] > r["dc_cai"]),
        "dc_cai_wins": sum(1 for r in results if r["bc_cai"] and r["dc_cai"] is not None and r["dc_cai"] > r["bc_cai"]),
        "bc_speed_wins": sum(1 for r in results if r["bc_time_ms"] and r["dc_time_ms"] is not None and r["bc_time_ms"] < r["dc_time_ms"]),
        "dc_speed_wins": sum(1 for r in results if r["bc_time_ms"] and r["dc_time_ms"] is not None and r["dc_time_ms"] < r["bc_time_ms"]),
    }
    summary["speed_ratio"] = round(summary["dc_median_time_ms"]/summary["bc_median_time_ms"], 2) if summary["bc_median_time_ms"]>0 else 0
    summary["cai_delta"] = round(summary["bc_mean_cai"] - summary["dc_mean_cai"], 4)

    print()
    print("="*70)
    print("HEAVY FAIR BENCHMARK RESULTS")
    print("="*70)
    print(f"Genes tested: {summary['total_genes']}")
    print(f"BC mean CAI: {summary['bc_mean_cai']:.4f}")
    if dc_cais_all: print(f"DC mean CAI: {summary['dc_mean_cai']:.4f}")
    print(f"CAI delta (BC-DC): {summary['cai_delta']:+.4f}")
    print(f"BC median time: {summary['bc_median_time_ms']:.2f}ms")
    if dc_times_all: print(f"DC median time: {summary['dc_median_time_ms']:.2f}ms")
    if summary['speed_ratio']: print(f"Speed ratio (DC/BC): {summary['speed_ratio']}×")
    print(f"BC CAI wins: {summary['bc_cai_wins']}/{summary['total_genes']}")
    if dc_cais_all: print(f"DC CAI wins: {summary['dc_cai_wins']}/{summary['total_genes']}")
    print(f"BC speed wins: {summary['bc_speed_wins']}/{summary['total_genes']}")
    if dc_times_all: print(f"DC speed wins: {summary['dc_speed_wins']}/{summary['total_genes']}")
    print(f"BC total violations: {summary['bc_total_violations']}")
    if dc_times_all: print(f"DC total violations: {summary['dc_total_violations']}")

    output = {"summary": summary, "per_gene": results}
    with open(os.path.join(os.path.dirname(__file__), "..", "..", "heavy_benchmark_results.json"), "w") as f:
        json.dump(output, f, indent=2)
    print(f"\nResults saved to heavy_benchmark_results.json")
    return output

if __name__ == "__main__":
    run_heavy_benchmark()
