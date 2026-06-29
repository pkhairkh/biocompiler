#!/usr/bin/env python3
"""Phase 3 incremental runner — resumes from existing CSV."""
import sys, time, gc, csv, json
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent.parent.parent / "src"))

OUTPUT = Path("/home/z/my-project/download/benchmark_all_from_scratch")
CSV_PATH = OUTPUT / "phase3_full_pipeline.csv"

from biocompiler.benchmarking.gene_sets import *
from biocompiler.shared.constants import AA_TO_CODONS
from biocompiler.expression.translation import compute_cai
from biocompiler.expression.tai import compute_tai
from biocompiler.organisms import CODON_ADAPTIVENESS_TABLES, resolve_organism
from biocompiler.optimizer.pipeline_core import optimize_sequence

# Collect ALL proteins
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

PIPE_ORGS = [('Homo_sapiens','Human','eukaryote'), ('Escherichia_coli','E.coli','prokaryote'),
             ('Saccharomyces_cerevisiae','Yeast','eukaryote'), ('CHO_K1','CHO-K1','eukaryote'),
             ('Komagataella_phaffii','Pichia','eukaryote')]

# Load existing
existing = []
done_keys = set()
if CSV_PATH.exists():
    with open(CSV_PATH, 'r', newline='') as f:
        existing = list(csv.DictReader(f))
        for r in existing:
            done_keys.add((r['protein_name'], r['organism']))

print(f"Phase 3: {len(proteins)} proteins × 5 organisms, {len(done_keys)} already done", flush=True)

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

results = list(existing)
new_count = 0

for ip, prot in enumerate(proteins):
    any_new = False
    for org_canonical, org_display, org_domain in PIPE_ORGS:
        key = (prot['name'], org_display)
        if key in done_keys:
            continue

        protein_seq = prot['protein']
        row = {'protein_name': prot['name'], 'protein_length': prot['length_aa'],
               'protein_category': prot['category'], 'organism': org_display, 'organism_domain': org_domain}

        # Naive
        try:
            nseq = naive_cai_seq(protein_seq, org_canonical)
            row['naive_cai_ceiling'] = round(compute_cai(nseq, organism=org_canonical), 4)
            row['naive_gc'] = round(compute_gc(nseq)*100, 1)
            row['naive_cpg'] = count_cpg(nseq)
            try: row['naive_tai'] = round(compute_tai(nseq, organism=org_canonical), 4)
            except: row['naive_tai'] = 0.0
        except: row.update({'naive_cai_ceiling':0.0, 'naive_gc':0.0, 'naive_cpg':0, 'naive_tai':0.0})

        # Pipeline
        t0 = time.perf_counter()
        try:
            bcr = optimize_sequence(target_protein=protein_seq, organism=org_canonical,
                gc_lo=0.30, gc_hi=0.70, strategy='hybrid',
                optimize_mrna_stability=False, include_utr=False,
                consider_codon_pair_bias=False, track_provenance=False, strict_mode=False)
            row['pipeline_cai'] = round(bcr.cai, 4)
            row['pipeline_gc'] = round(compute_gc(bcr.sequence)*100, 1)
            row['pipeline_cpg'] = count_cpg(bcr.sequence)
            try: row['pipeline_tai'] = round(compute_tai(bcr.sequence, organism=org_canonical), 4)
            except: row['pipeline_tai'] = 0.0
            row['cai_cost'] = round(row['naive_cai_ceiling'] - row['pipeline_cai'], 4) if row['naive_cai_ceiling'] > 0 else 0.0
            row['cpg_reduction'] = row['naive_cpg'] - row['pipeline_cpg']
            row['failed_predicates'] = ';'.join(bcr.failed_predicates) if bcr.failed_predicates else ''
            row['n_failed_predicates'] = len(bcr.failed_predicates) if bcr.failed_predicates else 0
            row['success'] = True; row['blocked'] = False; row['error'] = ''
        except Exception as e:
            err = str(e)
            if 'BIOSECURITY' in err.upper():
                row.update({'pipeline_cai':0.0,'pipeline_gc':0.0,'pipeline_cpg':0,'pipeline_tai':0.0,
                    'cai_cost':0.0,'cpg_reduction':0,'failed_predicates':'','n_failed_predicates':0,
                    'success':True,'blocked':True,'error':f'BLOCKED: {err[:60]}'})
            else:
                row.update({'pipeline_cai':0.0,'pipeline_gc':0.0,'pipeline_cpg':0,'pipeline_tai':0.0,
                    'cai_cost':0.0,'cpg_reduction':0,'failed_predicates':'','n_failed_predicates':0,
                    'success':False,'blocked':False,'error':err[:100]})
        row['time_s'] = round(time.perf_counter() - t0, 3)

        results.append(row)
        new_count += 1
        any_new = True

    if any_new:
        # Save incrementally
        if results:
            fieldnames = ['protein_name','protein_length','protein_category','organism','organism_domain',
                'naive_cai_ceiling','naive_tai','naive_gc','naive_cpg',
                'pipeline_cai','pipeline_tai','pipeline_gc','pipeline_cpg',
                'cai_cost','cpg_reduction','n_failed_predicates','failed_predicates','time_s','success','blocked','error']
            with open(CSV_PATH, 'w', newline='') as f:
                w = csv.DictWriter(f, fieldnames=fieldnames, extrasaction='ignore')
                w.writeheader(); w.writerows(results)
        print(f"  [{ip+1}/{len(proteins)}] {prot['name']} ({prot['length_aa']}aa) done", flush=True)
    gc.collect()

print(f"\nPhase 3 COMPLETE: {len(results)} total results ({new_count} new)", flush=True)
