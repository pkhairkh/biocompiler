#!/usr/bin/env python3
"""Phase 1: Head-to-head benchmark — runs incrementally, saves after each protein."""
import sys, time, gc, csv, json, os
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent.parent.parent / "src"))

OUTPUT = Path("/home/z/my-project/download/benchmark_all_from_scratch")
OUTPUT.mkdir(parents=True, exist_ok=True)
LOG = OUTPUT / "phase1.log"

def log(msg):
    with open(LOG, "a") as f:
        f.write(f"{time.strftime('%H:%M:%S')} {msg}\n")
    print(msg, flush=True)

from biocompiler.benchmarking.gene_sets import *
from biocompiler.shared.constants import AA_TO_CODONS
from biocompiler.expression.translation import compute_cai
from biocompiler.expression.tai import compute_tai
from biocompiler.organisms import CODON_ADAPTIVENESS_TABLES, resolve_organism
from biocompiler.benchmarking.dnachisel_adapter import DNAchiselAdapter
from biocompiler.optimizer.pipeline_core import optimize_sequence

# Collect ALL proteins, sorted by length
proteins = []
seen = set()
def _add(name, seq, source, cat):
    key = f'{name}_{source}'
    if key in seen: return
    seen.add(key)
    if any(aa not in 'ACDEFGHIKLMNPQRSTVWY*' for aa in seq): return
    proteins.append({'name': name, 'protein': seq, 'source': source, 'category': cat, 'length_aa': len(seq)})
for n, d in HUMAN_THERAPEUTIC_GENES.items(): _add(n, d['protein_sequence'], 'HT', 'therapeutic')
for n, d in VACCINE_ANTIGEN_GENES.items(): _add(n, d['protein_sequence'], 'VA', 'vaccine')
for n, d in E_COLI_EXTENDED.items(): _add(n, d['protein'], 'EC', d.get('category', 'ecoli'))
for n, d in HUMAN_SIGNALING.items(): _add(n, d['protein'], 'HS', d.get('category', 'signaling'))
for n, d in YEAST_INDUSTRIAL.items(): _add(n, d['protein'], 'YI', d.get('category', 'industrial'))
for n, d in MOUSE_MODEL.items(): _add(n, d['protein'], 'MM', d.get('category', 'model_organism'))
for n, (seq, desc) in BENCHMARK_GENES.get('standard', {}).items(): _add(n, seq, 'BS', 'standard')
proteins.sort(key=lambda p: p['length_aa'])

# Load existing results for resume
csv_path = OUTPUT / "phase1_h2h.csv"
existing = []
done_keys = set()
if csv_path.exists():
    with open(csv_path, 'r', newline='') as f:
        reader = csv.DictReader(f)
        existing = list(reader)
        for r in existing:
            done_keys.add((r['protein_name'], r['organism']))
    log(f"Resuming: {len(done_keys)} already done")

def compute_gc(s): return sum(1 for b in s.upper() if b in 'GC') / len(s) if s else 0.0
def count_cpg(s): return s.upper().count('CG') if s else 0
def naive_cai_seq(protein, organism):
    resolved = resolve_organism(organism, strict=False)
    adaptiveness = CODON_ADAPTIVENESS_TABLES.get(resolved, CODON_ADAPTIVENESS_TABLES.get('Homo_sapiens'))
    codons = []
    for aa in protein:
        if aa == 'M': codons.append('ATG'); continue
        if aa == '*': codons.append('TAA'); continue
        cands = AA_TO_CODONS.get(aa, [])
        if not cands: codons.append('NNN'); continue
        codons.append(max(cands, key=lambda c: adaptiveness.get(c, 0.0)))
    return ''.join(codons)

H2H_ORGS = [('Homo_sapiens','Human'), ('Escherichia_coli','E.coli'), ('Saccharomyces_cerevisiae','Yeast'), ('Mus_musculus','Mouse')]
adapter = DNAchiselAdapter()

results = list(existing)
new_count = 0

log(f"PHASE 1: {len(proteins)} proteins × 4 organisms ({len(done_keys)} already done)")

for ip, prot in enumerate(proteins):
    for org_canonical, org_display in H2H_ORGS:
        key = (prot['name'], org_display)
        if key in done_keys:
            continue

        row = {'protein_name': prot['name'], 'length': prot['length_aa'], 'organism': org_display, 'category': prot['category']}
        protein_seq = prot['protein']

        # Naive
        try:
            nseq = naive_cai_seq(protein_seq, org_canonical)
            row['naive_cai'] = round(compute_cai(nseq, organism=org_canonical), 4)
            row['naive_cpg'] = count_cpg(nseq)
            row['naive_gc'] = round(compute_gc(nseq)*100, 1)
            try: row['naive_tai'] = round(compute_tai(nseq, organism=org_canonical), 4)
            except: row['naive_tai'] = 0.0
        except: row.update({'naive_cai':0.0, 'naive_cpg':0, 'naive_gc':0.0, 'naive_tai':0.0})

        # DNAchisel
        t0 = time.perf_counter()
        try:
            dcr = adapter.optimize(protein_seq, organism=org_canonical, constraints=[{'type':'gc_range','gc_lo':0.30,'gc_hi':0.70}])
            row['dc_cai'] = round(dcr.cai, 4)
            row['dc_cpg'] = count_cpg(dcr.sequence) if dcr.sequence else 0
            row['dc_gc'] = round(compute_gc(dcr.sequence)*100,1) if dcr.sequence else 0.0
            row['dc_success'] = dcr.success
            try: row['dc_tai'] = round(compute_tai(dcr.sequence, organism=org_canonical), 4) if dcr.sequence else 0.0
            except: row['dc_tai'] = 0.0
        except Exception as e:
            row.update({'dc_cai':0.0, 'dc_cpg':0, 'dc_gc':0.0, 'dc_success':False, 'dc_tai':0.0})
        row['dc_time'] = round(time.perf_counter()-t0, 2)

        # BioCompiler
        t0 = time.perf_counter()
        try:
            bcr = optimize_sequence(target_protein=protein_seq, organism=org_canonical,
                gc_lo=0.30, gc_hi=0.70, strategy='hybrid',
                optimize_mrna_stability=False, include_utr=False,
                consider_codon_pair_bias=False, track_provenance=False, strict_mode=False)
            row['bc_cai'] = round(bcr.cai, 4)
            row['bc_cpg'] = count_cpg(bcr.sequence)
            row['bc_gc'] = round(compute_gc(bcr.sequence)*100, 1)
            row['bc_success'] = True; row['bc_blocked'] = False
            try: row['bc_tai'] = round(compute_tai(bcr.sequence, organism=org_canonical), 4)
            except: row['bc_tai'] = 0.0
            row['failed_preds'] = ';'.join(bcr.failed_predicates) if bcr.failed_predicates else ''
            row['n_failed'] = len(bcr.failed_predicates) if bcr.failed_predicates else 0
        except Exception as e:
            err = str(e)
            if 'BIOSECURITY' in err.upper():
                row.update({'bc_cai':0.0,'bc_cpg':0,'bc_gc':0.0,'bc_success':True,'bc_blocked':True,'bc_tai':0.0,'failed_preds':'','n_failed':0})
            else:
                row.update({'bc_cai':0.0,'bc_cpg':0,'bc_gc':0.0,'bc_success':False,'bc_blocked':False,'bc_tai':0.0,'failed_preds':'','n_failed':0,'bc_err':err[:80]})
        row['bc_time'] = round(time.perf_counter()-t0, 2)
        row['cai_delta'] = round(row.get('bc_cai',0) - row.get('dc_cai',0), 4)
        row['cpg_delta'] = row.get('bc_cpg',0) - row.get('dc_cpg',0)
        results.append(row)
        new_count += 1

    log(f"  [{ip+1}/{len(proteins)}] {prot['name']} ({prot['length_aa']}aa)")

    # Save every protein
    if results:
        with open(csv_path, 'w', newline='') as f:
            w = csv.DictWriter(f, fieldnames=list(results[0].keys()))
            w.writeheader(); w.writerows(results)
    gc.collect()

log(f"PHASE 1 COMPLETE: {len(results)} total results ({new_count} new)")
