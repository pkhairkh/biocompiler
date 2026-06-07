#!/usr/bin/env python3
"""
Generate benchmark visualization charts from e2e_benchmark_results.json
"""

import json
import statistics
from pathlib import Path
from collections import defaultdict

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.font_manager as fm
import numpy as np

# Font setup
fm.fontManager.addfont('/usr/share/fonts/truetype/chinese/SarasaMonoSC-Regular.ttf')
fm.fontManager.addfont('/usr/share/fonts/truetype/dejavu/DejaVuSans.ttf')
plt.rcParams['font.sans-serif'] = ['DejaVu Sans', 'Sarasa Mono SC']
plt.rcParams['axes.unicode_minus'] = False

# Load data
with open("/home/z/my-project/biocompiler/benchmark_results/e2e_benchmark_results.json") as f:
    data = json.load(f)

output_dir = Path("/home/z/my-project/biocompiler/benchmark_results")

# Color palette
COLORS = {
    'naive': '#94A3B8',       # gray
    'biocompiler': '#CC3311', # red (alert)
    'dnachisel': '#0077BB',   # blue
    'dnachisel+re': '#33BBEE', # light blue
}

# ═══════════════════════════════════════════════════════════
# Chart 1: Head-to-Head CAI Comparison (E. coli genes)
# ═══════════════════════════════════════════════════════════

ecoli_genes = [r for r in data if r['organism'] == 'Escherichia_coli' and r['tool'] in ('biocompiler', 'dnachisel', 'naive') and not r['error']]

gene_tools = defaultdict(dict)
for r in ecoli_genes:
    gene_tools[r['gene']][r['tool']] = r['cais'].get('Escherichia_coli', 0)

genes = sorted(set(k for k, v in gene_tools.items() if 'biocompiler' in v and 'dnachisel' in v))
bc_cais = [gene_tools[g].get('biocompiler', 0) for g in genes]
dc_cais = [gene_tools[g].get('dnachisel', 0) for g in genes]
nv_cais = [gene_tools[g].get('naive', 0) for g in genes]

fig, ax = plt.subplots(figsize=(16, 7))
x = np.arange(len(genes))
width = 0.25

bars1 = ax.bar(x - width, nv_cais, width, label='Naive (most-freq codon)', color=COLORS['naive'], edgecolor='white', linewidth=0.5)
bars2 = ax.bar(x, bc_cais, width, label='BioCompiler', color=COLORS['biocompiler'], edgecolor='white', linewidth=0.5)
bars3 = ax.bar(x + width, dc_cais, width, label='DNAchisel', color=COLORS['dnachisel'], edgecolor='white', linewidth=0.5)

ax.set_xlabel('Gene', fontsize=13)
ax.set_ylabel('CAI (Codon Adaptation Index)', fontsize=13)
ax.set_title('E. coli Optimization: BioCompiler vs DNAchisel vs Naive Baseline', fontsize=16, fontweight='bold')
ax.set_xticks(x)
ax.set_xticklabels(genes, rotation=45, ha='right', fontsize=10)
ax.set_ylim(0, 1.1)
ax.legend(loc='upper right', fontsize=11, framealpha=0.9)
ax.axhline(y=0.8, color='#EE7733', linestyle='--', alpha=0.5, linewidth=1, label='CAI=0.8 target')
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
ax.grid(axis='y', alpha=0.1)

plt.tight_layout()
plt.savefig(output_dir / 'chart1_ecoli_cai_comparison.png', dpi=200, bbox_inches='tight')
plt.close()
print(f"Chart 1 saved: {output_dir / 'chart1_ecoli_cai_comparison.png'}")


# ═══════════════════════════════════════════════════════════
# Chart 2: CAI Delta (BioCompiler - DNAchisel) per gene
# ═══════════════════════════════════════════════════════════

all_genes_h2h = defaultdict(dict)
for r in data:
    if r['tool'] in ('biocompiler', 'dnachisel') and not r['error']:
        all_genes_h2h[r['gene']][r['tool']] = r['cais'].get(r['organism'], 0)

genes_all = sorted(all_genes_h2h.keys())
deltas = []
for g in genes_all:
    bc = all_genes_h2h[g].get('biocompiler', 0)
    dc_val = all_genes_h2h[g].get('dnachisel', 0)
    deltas.append(bc - dc_val)

fig, ax = plt.subplots(figsize=(16, 7))
colors_bar = ['#CC3311' if d < 0 else '#009988' for d in deltas]
bars = ax.bar(range(len(genes_all)), deltas, color=colors_bar, edgecolor='white', linewidth=0.5, width=0.8)

ax.set_xlabel('Gene', fontsize=13)
ax.set_ylabel('CAI Delta (BioCompiler - DNAchisel)', fontsize=13)
ax.set_title('CAI Gap: BioCompiler vs DNAchisel (Negative = BioCompiler Underperforms)', fontsize=15, fontweight='bold')
ax.set_xticks(range(len(genes_all)))
ax.set_xticklabels(genes_all, rotation=55, ha='right', fontsize=9)
ax.axhline(y=0, color='#37352F', linewidth=1.5)
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
ax.grid(axis='y', alpha=0.1)

# Add mean delta line
mean_delta = statistics.mean(deltas)
ax.axhline(y=mean_delta, color='#CC3311', linestyle='--', linewidth=1.5, alpha=0.7)
ax.text(len(genes_all)-1, mean_delta + 0.02, f'Mean = {mean_delta:.3f}', fontsize=11, color='#CC3311', ha='right')

plt.tight_layout()
plt.savefig(output_dir / 'chart2_cai_delta.png', dpi=200, bbox_inches='tight')
plt.close()
print(f"Chart 2 saved: {output_dir / 'chart2_cai_delta.png'}")


# ═══════════════════════════════════════════════════════════
# Chart 3: Speed Comparison (log scale)
# ═══════════════════════════════════════════════════════════

tools_speed = defaultdict(list)
for r in data:
    if not r['error'] and r['runtime_s'] > 0:
        tools_speed[r['tool']].append(r['runtime_s'] * 1000)  # ms

fig, ax = plt.subplots(figsize=(10, 6))
tool_names = ['naive', 'biocompiler', 'dnachisel', 'dnachisel+re']
tool_labels = ['Naive', 'BioCompiler', 'DNAchisel', 'DNAchisel+RE']

bp_data = [tools_speed.get(t, [0]) for t in tool_names]
bp = ax.boxplot(bp_data, labels=tool_labels, patch_artist=True, widths=0.5,
                medianprops={'color': '#37352F', 'linewidth': 2})

for patch, color_key in zip(bp['boxes'], tool_names):
    patch.set_facecolor(COLORS.get(color_key, '#E5E7EB'))
    patch.set_alpha(0.7)

ax.set_ylabel('Runtime (ms, log scale)', fontsize=13)
ax.set_title('Optimization Speed: BioCompiler vs DNAchisel', fontsize=16, fontweight='bold')
ax.set_yscale('log')
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
ax.grid(axis='y', alpha=0.1)

# Add mean annotations
for i, t in enumerate(tool_names):
    vals = tools_speed.get(t, [0])
    if vals and vals[0] > 0:
        ax.text(i+1, statistics.mean(vals)*1.5, f'mean={statistics.mean(vals):.1f}ms',
                ha='center', fontsize=9, color='#37352F')

plt.tight_layout()
plt.savefig(output_dir / 'chart3_speed_comparison.png', dpi=200, bbox_inches='tight')
plt.close()
print(f"Chart 3 saved: {output_dir / 'chart3_speed_comparison.png'}")


# ═══════════════════════════════════════════════════════════
# Chart 4: Constraint Violations in BioCompiler
# ═══════════════════════════════════════════════════════════

bc_results = [r for r in data if r['tool'] == 'biocompiler' and not r['error']]
violation_counts = defaultdict(int)
for r in bc_results:
    for v in r.get('constraint_violations', []):
        violation_counts[v] += 1

violations = sorted(violation_counts.items(), key=lambda x: -x[1])
v_names = [v[0] for v in violations]
v_counts = [v[1] for v in violations]
v_pcts = [c/len(bc_results)*100 for c in v_counts]

fig, ax = plt.subplots(figsize=(10, 6))
bars = ax.barh(range(len(v_names)), v_pcts, color='#CC3311', alpha=0.75, edgecolor='white', height=0.6)

ax.set_yticks(range(len(v_names)))
ax.set_yticklabels(v_names, fontsize=11)
ax.set_xlabel('% of Genes with Violation', fontsize=13)
ax.set_title('BioCompiler: Constraint Violations After Optimization', fontsize=15, fontweight='bold')
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
ax.grid(axis='x', alpha=0.1)

for i, (pct, count) in enumerate(zip(v_pcts, v_counts)):
    ax.text(pct + 1, i, f'{pct:.0f}% ({count}/{len(bc_results)})', va='center', fontsize=10, color='#37352F')

ax.set_xlim(0, max(v_pcts) * 1.25)
plt.tight_layout()
plt.savefig(output_dir / 'chart4_violations.png', dpi=200, bbox_inches='tight')
plt.close()
print(f"Chart 4 saved: {output_dir / 'chart4_violations.png'}")


# ═══════════════════════════════════════════════════════════
# Chart 5: Cross-organism EGFP CAI
# ═══════════════════════════════════════════════════════════

egfp_data = [r for r in data if r['gene'] == 'EGFP']
org_tools = defaultdict(dict)
for r in egfp_data:
    if not r['error']:
        org_tools[r['organism']][r['tool']] = r['cais'].get(r['organism'], 0)

organisms = ['Escherichia_coli', 'Saccharomyces_cerevisiae', 'Homo_sapiens', 'Mus_musculus', 'CHO_K1']
org_labels = ['E. coli', 'S. cerevisiae', 'H. sapiens', 'M. musculus', 'CHO-K1']

fig, ax = plt.subplots(figsize=(12, 6))
x = np.arange(len(organisms))
width = 0.25

nv_vals = [org_tools[o].get('naive', 0) for o in organisms]
bc_vals = [org_tools[o].get('biocompiler', 0) for o in organisms]
dc_vals = [org_tools[o].get('dnachisel', 0) for o in organisms]

ax.bar(x - width, nv_vals, width, label='Naive', color=COLORS['naive'], edgecolor='white')
ax.bar(x, bc_vals, width, label='BioCompiler', color=COLORS['biocompiler'], edgecolor='white')
ax.bar(x + width, dc_vals, width, label='DNAchisel', color=COLORS['dnachisel'], edgecolor='white')

ax.set_xticks(x)
ax.set_xticklabels(org_labels, fontsize=12)
ax.set_ylabel('CAI', fontsize=13)
ax.set_title('Cross-Organism EGFP Optimization: CAI Comparison', fontsize=15, fontweight='bold')
ax.set_ylim(0, 1.15)
ax.legend(loc='best', fontsize=11)
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
ax.grid(axis='y', alpha=0.1)

plt.tight_layout()
plt.savefig(output_dir / 'chart5_cross_organism.png', dpi=200, bbox_inches='tight')
plt.close()
print(f"Chart 5 saved: {output_dir / 'chart5_cross_organism.png'}")


# ═══════════════════════════════════════════════════════════
# Chart 6: GC Content Distribution
# ═══════════════════════════════════════════════════════════

fig, axes = plt.subplots(1, 2, figsize=(14, 6))

# E. coli GC
ecoli_tools = defaultdict(list)
for r in data:
    if r['organism'] == 'Escherichia_coli' and not r['error'] and r['tool'] in ('biocompiler', 'dnachisel'):
        ecoli_tools[r['tool']].append(r['gc_content'])

for tool, vals in ecoli_tools.items():
    axes[0].hist(vals, bins=10, alpha=0.6, label=tool, color=COLORS.get(tool, '#999'), edgecolor='white')

axes[0].set_xlabel('GC Content', fontsize=12)
axes[0].set_ylabel('Count', fontsize=12)
axes[0].set_title('E. coli: GC Distribution', fontsize=14, fontweight='bold')
axes[0].legend(fontsize=11)
axes[0].axvspan(0.30, 0.70, alpha=0.08, color='green', label='Target GC range')
axes[0].spines['top'].set_visible(False)
axes[0].spines['right'].set_visible(False)

# Human GC
human_tools = defaultdict(list)
for r in data:
    if r['organism'] == 'Homo_sapiens' and not r['error'] and r['tool'] in ('biocompiler', 'dnachisel'):
        human_tools[r['tool']].append(r['gc_content'])

for tool, vals in human_tools.items():
    axes[1].hist(vals, bins=8, alpha=0.6, label=tool, color=COLORS.get(tool, '#999'), edgecolor='white')

axes[1].set_xlabel('GC Content', fontsize=12)
axes[1].set_ylabel('Count', fontsize=12)
axes[1].set_title('Human: GC Distribution', fontsize=14, fontweight='bold')
axes[1].legend(fontsize=11)
axes[1].axvspan(0.30, 0.70, alpha=0.08, color='green')
axes[1].spines['top'].set_visible(False)
axes[1].spines['right'].set_visible(False)

plt.tight_layout()
plt.savefig(output_dir / 'chart6_gc_distribution.png', dpi=200, bbox_inches='tight')
plt.close()
print(f"Chart 6 saved: {output_dir / 'chart6_gc_distribution.png'}")


# ═══════════════════════════════════════════════════════════
# Chart 7: Summary Dashboard
# ═══════════════════════════════════════════════════════════

fig, axes = plt.subplots(2, 2, figsize=(14, 10))

# 7a: Mean CAI by tool and domain
prokaryote_genes = [r for r in data if r['organism'] == 'Escherichia_coli' and not r['error']]
eukaryote_genes = [r for r in data if r['organism'] in ('Homo_sapiens', 'Mus_musculus', 'CHO_K1', 'Saccharomyces_cerevisiae') and not r['error']]

tools_summary = ['naive', 'biocompiler', 'dnachisel']
prok_cai = [statistics.mean([r['cais'].get(r['organism'], 0) for r in prokaryote_genes if r['tool'] == t]) for t in tools_summary]
euk_cai = [statistics.mean([r['cais'].get(r['organism'], 0) for r in eukaryote_genes if r['tool'] == t]) for t in tools_summary]

x = np.arange(len(tools_summary))
width = 0.3
axes[0,0].bar(x - width/2, prok_cai, width, label='Prokaryote (E. coli)', color='#0077BB', alpha=0.8)
axes[0,0].bar(x + width/2, euk_cai, width, label='Eukaryote', color='#009988', alpha=0.8)
axes[0,0].set_xticks(x)
axes[0,0].set_xticklabels(['Naive', 'BioCompiler', 'DNAchisel'], fontsize=10)
axes[0,0].set_ylabel('Mean CAI')
axes[0,0].set_title('Mean CAI by Domain', fontweight='bold')
axes[0,0].set_ylim(0, 1.1)
axes[0,0].legend(fontsize=9)
axes[0,0].spines['top'].set_visible(False)
axes[0,0].spines['right'].set_visible(False)

# 7b: Speed comparison
tools_speed_bc = [r['runtime_s']*1000 for r in data if r['tool'] == 'biocompiler' and not r['error']]
tools_speed_dc = [r['runtime_s']*1000 for r in data if r['tool'] == 'dnachisel' and not r['error']]

axes[0,1].boxplot([tools_speed_bc, tools_speed_dc], labels=['BioCompiler', 'DNAchisel'],
                  patch_artist=True, widths=0.5)
axes[0,1].set_ylabel('Runtime (ms)')
axes[0,1].set_title('Speed: BioCompiler vs DNAchisel', fontweight='bold')
axes[0,1].spines['top'].set_visible(False)
axes[0,1].spines['right'].set_visible(False)

# Add speedup factor
bc_mean = statistics.mean(tools_speed_bc)
dc_mean = statistics.mean(tools_speed_dc)
speedup = bc_mean / dc_mean if dc_mean > 0 else 0
axes[0,1].text(0.5, 0.95, f'{speedup:.0f}x slower', transform=axes[0,1].transAxes,
              ha='center', va='top', fontsize=14, fontweight='bold', color='#CC3311')

# 7c: Constraint violations pie
viol_items = sorted(violation_counts.items(), key=lambda x: -x[1])[:6]
v_labels = [v[0] for v in viol_items]
v_vals = [v[1] for v in viol_items]
colors_pie = ['#CC3311', '#EE7733', '#0077BB', '#33BBEE', '#009988', '#999999']
axes[1,0].pie(v_vals, labels=v_labels, colors=colors_pie[:len(v_vals)], autopct='%1.0f%%',
             startangle=90, textprops={'fontsize': 9})
axes[1,0].set_title('Constraint Violations\n(BioCompiler)', fontweight='bold', fontsize=12)

# 7d: CAI distribution violin-like
bc_all_cai = [r['cais'].get(r['organism'], 0) for r in data if r['tool'] == 'biocompiler' and not r['error']]
dc_all_cai = [r['cais'].get(r['organism'], 0) for r in data if r['tool'] == 'dnachisel' and not r['error']]

axes[1,1].hist(bc_all_cai, bins=15, alpha=0.6, label=f'BioCompiler (mean={statistics.mean(bc_all_cai):.3f})',
              color=COLORS['biocomplier'] if 'biocomplier' in COLORS else COLORS['biocompiler'], edgecolor='white')
axes[1,1].hist(dc_all_cai, bins=15, alpha=0.6, label=f'DNAchisel (mean={statistics.mean(dc_all_cai):.3f})',
              color=COLORS['dnachisel'], edgecolor='white')
axes[1,1].set_xlabel('CAI')
axes[1,1].set_ylabel('Count')
axes[1,1].set_title('CAI Distribution', fontweight='bold')
axes[1,1].legend(fontsize=9)
axes[1,1].spines['top'].set_visible(False)
axes[1,1].spines['right'].set_visible(False)

plt.suptitle('BioCompiler E2E Benchmark Dashboard', fontsize=18, fontweight='bold', y=1.02)
plt.tight_layout()
plt.savefig(output_dir / 'chart7_dashboard.png', dpi=200, bbox_inches='tight')
plt.close()
print(f"Chart 7 saved: {output_dir / 'chart7_dashboard.png'}")

print("\nAll charts generated!")
