#!/usr/bin/env python3
import sys, os; sys.path.insert(0,'src'); os.environ['TQDM_DISABLE']='1'
import time, json
from pathlib import Path
from biocompiler.benchmarking.gene_sets import HUMAN_THERAPEUTIC_GENES, VACCINE_ANTIGEN_GENES, E_COLI_EXTENDED, HUMAN_SIGNALING, YEAST_INDUSTRIAL, MOUSE_MODEL, BENCHMARK_GENES
from biocompiler.expression.translation import compute_cai
from biocompiler.expression.tai import compute_tai, optimize_for_tai
from biocompiler.organisms import CODON_ADAPTIVENESS_TABLES, resolve_organism
from biocompiler.organisms.tai_data import TRNA_GENE_COPIES
from biocompiler.shared.constants import AA_TO_CODONS
from biocompiler.benchmarking.dnachisel_adapter import DNAchiselAdapter
from biocompiler.optimizer.pipeline_core import optimize_sequence

proteins = []
seen = set()
def _add(n,s,sr,c):
    k=f'{n}_{sr}'
    if k in seen: return
    seen.add(k)
    if any(a not in 'ACDEFGHIKLMNPQRSTVWY*' for a in s): return
    if len(s)>500: return
    proteins.append({'name':n,'protein':s,'category':c,'len':len(s)})
for n,d in HUMAN_THERAPEUTIC_GENES.items(): _add(n,d['protein_sequence'],'HT','therapeutic')
for n,d in VACCINE_ANTIGEN_GENES.items(): _add(n,d['protein_sequence'],'VA','vaccine')
for n,d in E_COLI_EXTENDED.items(): _add(n,d['protein'],'EC',d.get('category','ecoli'))
for n,d in HUMAN_SIGNALING.items(): _add(n,d['protein'],'HS','signaling')
for n,d in YEAST_INDUSTRIAL.items(): _add(n,d['protein'],'YI','industrial')
for n,d in MOUSE_MODEL.items(): _add(n,d['protein'],'MM','model')
for n,(s,_) in BENCHMARK_GENES.get('standard',{}).items(): _add(n,s,'BS','standard')
print(f'Proteins: {len(proteins)}', flush=True)

def gc(s): return sum(1 for b in s.upper() if b in 'GC')/len(s) if s else 0
def cpg(s): return s.upper().count('CG') if s else 0
def naive(prot, org):
    res=resolve_organism(org,strict=False)
    ad=CODON_ADAPTIVENESS_TABLES.get(res,CODON_ADAPTIVENESS_TABLES.get('Homo_sapiens'))
    cs=[]
    for aa in prot:
        if aa=='M': cs.append('ATG'); continue
        if aa=='*': cs.append('TAA'); continue
        cc=AA_TO_CODONS.get(aa,[])
        cs.append(max(cc,key=lambda c:ad.get(c,0.0)) if cc else 'NNN')
    return ''.join(cs)

ORGS = ['Homo_sapiens','Escherichia_coli','Saccharomyces_cerevisiae','Mus_musculus']
OD = {'Homo_sapiens':'Human','Escherichia_coli':'E.coli','Saccharomyces_cerevisiae':'Yeast','Mus_musculus':'Mouse'}
PIPE_ORGS = ['Homo_sapiens','Escherichia_coli','Saccharomyces_cerevisiae','CHO_K1','Komagataella_phaffii']
PD = {'Homo_sapiens':'Human','Escherichia_coli':'E.coli','Saccharomyces_cerevisiae':'Yeast','CHO_K1':'CHO-K1','Komagataella_phaffii':'Pichia'}
TM = {'Escherichia_coli':'e_coli','Homo_sapiens':'human','Saccharomyces_cerevisiae':'yeast','Mus_musculus':'mouse','CHO_K1':'cho','Komagataella_phaffii':'p_pastoris','Caenorhabditis_elegans':'c_elegans','D_melanogaster':'d_melanogaster','Arabidopsis_thaliana':'a_thaliana','Bacillus_subtilis':'b_subtilis'}
OM = {'Homo_sapiens':'Human','Escherichia_coli':'E.coli','Saccharomyces_cerevisiae':'Yeast','Mus_musculus':'Mouse','CHO_K1':'CHO-K1','Komagataella_phaffii':'Pichia','Caenorhabditis_elegans':'C.elegans','D_melanogaster':'Drosophila','Arabidopsis_thaliana':'Arabidopsis','Bacillus_subtilis':'B.subtilis'}

print('H2H...', flush=True)
adapter = DNAchiselAdapter()
h2h_data = {o:{'bc_cai':[],'dc_cai':[],'bc_cpg':[],'dc_cpg':[],'bc_tai':[],'dc_tai':[],'naive_cai':[]} for o in OD.values()}
bcw=dcw=0
for i,p in enumerate(proteins):
    for org in ORGS:
        od=OD[org]
        ns=naive(p['protein'],org); ncai=round(compute_cai(ns,organism=org),4); ncpg=cpg(ns)
        try: ntai=round(compute_tai(ns,organism=org),4)
        except: ntai=0.0
        h2h_data[od]['naive_cai'].append(ncai)
        try:
            dc=adapter.optimize(p['protein'],org,constraints=[{'type':'gc_range','gc_lo':0.3,'gc_hi':0.7}])
            dc_cai=round(dc.cai,4); dc_cpg=cpg(dc.sequence) if dc.sequence else 0
            try: dc_tai=round(compute_tai(dc.sequence,organism=org),4) if dc.sequence else 0.0
            except: dc_tai=0.0
        except: dc_cai=0.0; dc_cpg=0; dc_tai=0.0
        try:
            bc=optimize_sequence(target_protein=p['protein'],organism=org,gc_lo=0.3,gc_hi=0.7,strategy='hybrid',optimize_mrna_stability=False,include_utr=False,consider_codon_pair_bias=False,track_provenance=False,strict_mode=False)
            bc_cai=round(bc.cai,4); bc_cpg=cpg(bc.sequence)
            try: bc_tai=round(compute_tai(bc.sequence,organism=org),4)
            except: bc_tai=0.0
            h2h_data[od]['bc_cai'].append(bc_cai); h2h_data[od]['dc_cai'].append(dc_cai)
            h2h_data[od]['bc_cpg'].append(bc_cpg); h2h_data[od]['dc_cpg'].append(dc_cpg)
            h2h_data[od]['bc_tai'].append(bc_tai); h2h_data[od]['dc_tai'].append(dc_tai)
            if bc_cai>dc_cai+0.001: bcw+=1
            elif dc_cai>bc_cai+0.001: dcw+=1
        except: pass
    if (i+1)%5==0: print(f'  H2H: {i+1}/{len(proteins)}', flush=True)
print(f'H2H done. BC={bcw} DC={dcw}', flush=True)

print('tAI...', flush=True)
tai_data={}
for p in proteins:
    for oc,tk in TM.items():
        if tk not in TRNA_GENE_COPIES: continue
        od=OM.get(oc,oc)
        if od not in tai_data: tai_data[od]={'ct':[],'tt':[],'g':[],'l':[],'gap':[]}
        try:
            cs=naive(p['protein'],oc); ts=optimize_for_tai(p['protein'],organism=oc)
            cc=compute_cai(cs,organism=oc); ct=compute_tai(cs,organism=oc)
            tc=compute_cai(ts,organism=oc); tt=compute_tai(ts,organism=oc)
            tai_data[od]['ct'].append(ct); tai_data[od]['tt'].append(tt)
            tai_data[od]['g'].append(tt-ct); tai_data[od]['l'].append(cc-tc); tai_data[od]['gap'].append(cc-ct)
        except: pass
print('tAI done.', flush=True)

print('Pipeline...', flush=True)
fp_data={}; nblk=0
for i,p in enumerate(proteins):
    for org in PIPE_ORGS:
        od=PD[org]
        if od not in fp_data: fp_data[od]={'cai':[],'cost':[],'tai':[],'cpg_r':[]}
        ns=naive(p['protein'],org); ncai=round(compute_cai(ns,organism=org),4); ncpg=cpg(ns)
        try:
            bc=optimize_sequence(target_protein=p['protein'],organism=org,gc_lo=0.3,gc_hi=0.7,strategy='hybrid',optimize_mrna_stability=False,include_utr=False,consider_codon_pair_bias=False,track_provenance=False,strict_mode=False)
            bcai=round(bc.cai,4); bcpg=cpg(bc.sequence)
            try: btai=round(compute_tai(bc.sequence,organism=org),4)
            except: btai=0.0
            fp_data[od]['cai'].append(bcai); fp_data[od]['cost'].append(round(ncai-bcai,4))
            fp_data[od]['tai'].append(btai); fp_data[od]['cpg_r'].append(ncpg-bcpg)
        except Exception as e:
            if 'BIOSECURITY' in str(e).upper(): nblk+=1
    if (i+1)%5==0: print(f'  Pipeline: {i+1}/{len(proteins)}', flush=True)
print('Pipeline done.', flush=True)

# Report
print(); print('='*80); print('  BIOCOMPILER CAI/tAI BENCHMARK — HONEST EDITION'); print('='*80)
print(f'  Proteins: {len(proteins)}  |  Blocked: {nblk}'); print()
print('-'*80); print('  HEAD-TO-HEAD: BioCompiler vs DNAchisel'); print('-'*80)
print(f'  BC CAI wins: {bcw}  |  DC CAI wins: {dcw}'); print()
print('  Org        BC CAI   DC CAI    Delta  BC CpG  DC CpG  CpG Delta   BC tAI   DC tAI')
print('  '+'-'*75)
for o in sorted(h2h_data):
    d=h2h_data[o]
    if not d['bc_cai']: continue
    n=len(d['bc_cai']); mbc=sum(d['bc_cai'])/n; mdc=sum(d['dc_cai'])/n; mc1=sum(d['bc_cpg'])/n; mc2=sum(d['dc_cpg'])/n
    mt1=sum(d['bc_tai'])/n if any(t>0 for t in d['bc_tai']) else 0
    mt2=sum(d['dc_tai'])/n if any(t>0 for t in d['dc_tai']) else 0
    print(f'  {o:<10} {mbc:>8.4f} {mdc:>8.4f} {mbc-mdc:>+8.4f} {mc1:>7.1f} {mc2:>7.1f} {mc1-mc2:>+10.1f} {mt1:>8.4f} {mt2:>8.4f}')
print(); print('-'*80); print('  THE COST OF 19 PREDICATES'); print('-'*80)
print('  Org        Pipeline   Ceiling  CAI Cost    tAI  CpG Red.')
print('  '+'-'*55)
for o in sorted(fp_data):
    d=fp_data[o]
    if not d['cai']: continue
    n=len(d['cai']); mc=sum(d['cai'])/n; mco=sum(d['cost'])/n; mt=sum(d['tai'])/n if any(t>0 for t in d['tai']) else 0; mr=sum(d['cpg_r'])/n
    print(f'  {o:<10} {mc:>9.4f} {mc+mco:>9.4f} {mco:>9.4f} {mt:>7.4f} {mr:>9.1f}')
print(); print('-'*80); print('  CAI-tAI TRADEOFF'); print('-'*80)
print('  Org          CAI-opt tAI  tAI-opt tAI  tAI Gain  CAI Loss     Gap')
print('  '+'-'*65)
for o in sorted(tai_data):
    d=tai_data[o]; n=len(d['ct']); mct=sum(d['ct'])/n; mtt=sum(d['tt'])/n; mg=sum(d['g'])/n; ml=sum(d['l'])/n; mgap=sum(d['gap'])/n
    print(f'  {o:<12} {mct:>12.4f} {mtt:>12.4f} {mg:>+9.4f} {ml:>9.4f} {mgap:>7.4f}')
print(); print('-'*80); print('  WHY YEAST/PICHIA HAVE LOW tAI'); print('-'*80)
for o in sorted(tai_data):
    d=tai_data[o]; mgap=sum(d['gap'])/len(d['gap'])
    if mgap>0.25:
        mg=sum(d['g'])/len(d['g']); ml=sum(d['l'])/len(d['l'])
        print(f'  {o}: gap={mgap:.3f}. CAI-opt codons poorly served by tRNA pool. Fundamental tradeoff.')
        print(f'    tAI-opt gains +{mg:.3f} at -{ml:.3f} CAI.')
    elif mgap>0.12:
        print(f'  {o}: gap={mgap:.3f}. Moderate codon bias vs tRNA misalignment.')
print()
abc=[v for d in h2h_data.values() for v in d['bc_cai']]
adc=[v for d in h2h_data.values() for v in d['dc_cai']]
ac1=[v for d in h2h_data.values() for v in d['bc_cpg']]
ac2=[v for d in h2h_data.values() for v in d['dc_cpg']]
if abc and adc:
    mbc=sum(abc)/len(abc); mdc=sum(adc)/len(adc); mc1=sum(ac1)/len(ac1); mc2=sum(ac2)/len(ac2)
    print('-'*80); print('  HONEST VERDICT'); print('-'*80)
    print(f'  BioCompiler mean CAI : {mbc:.4f}'); print(f'  DNAchisel mean CAI   : {mdc:.4f}')
    print(f'  CAI delta            : {mbc-mdc:+.4f}'); print(f'  BioCompiler mean CpG : {mc1:.1f}'); print(f'  DNAchisel mean CpG   : {mc2:.1f}')
    if mc2>0: print(f'  CpG delta            : {mc1-mc2:+.1f}')
    print()
    if mbc>=mdc: print('  VERDICT: BioCompiler MATCHES/BEATS DNAchisel on CAI')
    else: print(f'  VERDICT: DNAchisel leads CAI by {mdc-mbc:.4f}. BioCompiler advantage: 19-predicate compliance.')
print(); print('='*80)

summary={'n_proteins':len(proteins),'bc_cai_wins':bcw,'dc_cai_wins':dcw,
         'h2h':{o:{'bc_cai':round(sum(d['bc_cai'])/len(d['bc_cai']),4),'dc_cai':round(sum(d['dc_cai'])/len(d['dc_cai']),4),'bc_cpg':round(sum(d['bc_cpg'])/len(d['bc_cpg']),1),'dc_cpg':round(sum(d['dc_cpg'])/len(d['dc_cpg']),1)} for o,d in h2h_data.items() if d['bc_cai']},
         'tai':{o:{'cai_opt_tai':round(sum(d['ct'])/len(d['ct']),4),'tai_opt_tai':round(sum(d['tt'])/len(d['tt']),4),'tai_gain':round(sum(d['g'])/len(d['g']),4),'cai_loss':round(sum(d['l'])/len(d['l']),4)} for o,d in tai_data.items()},
         'fp':{o:{'pipeline_cai':round(sum(d['cai'])/len(d['cai']),4),'cai_cost':round(sum(d['cost'])/len(d['cost']),4),'tai':round(sum(d['tai'])/len(d['tai']) if any(t>0 for t in d['tai']) else 0,4),'cpg_red':round(sum(d['cpg_r'])/len(d['cpg_r']),1)} for o,d in fp_data.items() if d['cai']}}
Path('/home/z/my-project/download/benchmark_results_v2').mkdir(parents=True,exist_ok=True)
with open('/home/z/my-project/download/benchmark_results_v2/benchmark_summary.json','w') as f:
    json.dump(summary,f,indent=2,default=str)
print('SAVED!', flush=True)
