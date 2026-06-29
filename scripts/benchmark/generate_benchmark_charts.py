#!/usr/bin/env python3
"""
Generate publication-quality charts from the multi-organism CAI/tAI benchmark results.
"""
import json
import os
import sys
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.font_manager as fm
import numpy as np

# ── Font setup ──
fm.fontManager.addfont('/usr/share/fonts/truetype/chinese/SarasaMonoSC-Regular.ttf')
fm.fontManager.addfont('/usr/share/fonts/truetype/dejavu/DejaVuSans.ttf')
plt.rcParams['font.sans-serif'] = ['Sarasa Mono SC', 'DejaVu Sans']
plt.rcParams['axes.unicode_minus'] = False

# ── Color palette ──
PALETTE = {
    'primary': '#4C6EF5',
    'secondary': '#7C3AED',
    'accent': '#3AAFA9',
    'warm': '#EE7733',
    'danger': '#CC3311',
    'neutral': '#64748B',
    'bg': '#F8FAFC',
    'text': '#243447',
    'grid': '#E2E8F0',
}

ORGANISM_COLORS = {
    'Human': '#4C6EF5',
    'Mouse': '#7C3AED',
    'CHO-K1': '#3AAFA9',
    'Drosophila': '#EE7733',
    'B. subtilis': '#CC3311',
    'Arabidopsis': '#009988',
    'E. coli': '#0077BB',
    'C. elegans': '#EE3377',
    'Pichia': '#33BBEE',
    'Yeast': '#DDCC77',
}

OUTPUT_DIR = Path('/home/z/my-project/download/benchmark_charts')
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

DATA_DIR = Path('/home/z/my-project/download/benchmark_results')


def load_data():
    """Load benchmark results from JSON."""
    with open(DATA_DIR / 'full_results.json') as f:
        return json.load(f)


def chart_1_tai_by_organism(data):
    """Bar chart: mean tAI by organism (CAI-optimized sequences)."""
    tai_results = data['tai_results']
    
    # Group by organism
    org_data = {}
    for r in tai_results:
        if 'error' in r:
            continue
        org = r['organism_display']
        if org not in org_data:
            org_data[org] = {'cai': [], 'tai': [], 'gap': []}
        org_data[org]['cai'].append(r['cai'])
        org_data[org]['tai'].append(r['tai'])
        org_data[org]['gap'].append(r['cai_tai_gap'])
    
    organisms = sorted(org_data.keys(), key=lambda x: np.mean(org_data[x]['tai']), reverse=True)
    mean_tais = [np.mean(org_data[o]['tai']) for o in organisms]
    mean_cais = [np.mean(org_data[o]['cai']) for o in organisms]
    mean_gaps = [np.mean(org_data[o]['gap']) for o in organisms]
    
    fig, ax = plt.subplots(figsize=(14, 7), facecolor=PALETTE['bg'])
    ax.set_facecolor(PALETTE['bg'])
    
    x = np.arange(len(organisms))
    width = 0.35
    
    bars_cai = ax.bar(x - width/2, mean_cais, width, label='Mean CAI',
                      color=PALETTE['primary'], alpha=0.85, edgecolor='white', linewidth=0.5)
    bars_tai = ax.bar(x + width/2, mean_tais, width, label='Mean tAI',
                      color=PALETTE['accent'], alpha=0.85, edgecolor='white', linewidth=0.5)
    
    # Add value labels
    for bar in bars_cai:
        height = bar.get_height()
        ax.text(bar.get_x() + bar.get_width()/2., height + 0.008,
                f'{height:.3f}', ha='center', va='bottom', fontsize=9, color=PALETTE['text'])
    for bar in bars_tai:
        height = bar.get_height()
        ax.text(bar.get_x() + bar.get_width()/2., height + 0.008,
                f'{height:.3f}', ha='center', va='bottom', fontsize=9, color=PALETTE['text'])
    
    ax.set_xlabel('Organism', fontsize=12, color=PALETTE['text'], labelpad=10)
    ax.set_ylabel('Index Value', fontsize=12, color=PALETTE['text'], labelpad=10)
    ax.set_title('CAI vs tAI Across 10 Organisms (59 Proteins, CAI-Optimized Sequences)',
                 fontsize=14, fontweight='bold', color=PALETTE['text'], pad=16)
    ax.set_xticks(x)
    ax.set_xticklabels(organisms, rotation=30, ha='right', fontsize=10)
    ax.set_ylim(0, 1.12)
    ax.legend(loc='upper right', fontsize=11, framealpha=0.9)
    ax.grid(axis='y', alpha=0.15, color=PALETTE['grid'])
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.spines['left'].set_color(PALETTE['grid'])
    ax.spines['bottom'].set_color(PALETTE['grid'])
    
    plt.tight_layout()
    fig.savefig(OUTPUT_DIR / 'chart1_cai_vs_tai_by_organism.png', dpi=300, bbox_inches='tight',
                facecolor=PALETTE['bg'])
    plt.close(fig)
    print(f"  Saved: chart1_cai_vs_tai_by_organism.png")


def chart_2_tai_gap_heatmap(data):
    """Heatmap: CAI-tAI gap per protein per organism."""
    tai_results = data['tai_results']
    
    # Get organisms and proteins
    organisms = []
    proteins = []
    gap_matrix = {}
    
    for r in tai_results:
        if 'error' in r:
            continue
        org = r['organism_display']
        prot = r['protein_name']
        if org not in organisms:
            organisms.append(org)
        if prot not in proteins:
            proteins.append(prot)
        gap_matrix[(prot, org)] = r['cai_tai_gap']
    
    # Sort organisms by mean gap
    org_mean_gap = {}
    for org in organisms:
        gaps = [gap_matrix.get((p, org), 0) for p in proteins if (p, org) in gap_matrix]
        org_mean_gap[org] = np.mean(gaps) if gaps else 0
    organisms = sorted(organisms, key=lambda x: org_mean_gap[x])
    
    # Select representative proteins (top 30 by variety)
    selected = proteins[:30]
    
    # Build matrix
    matrix = np.zeros((len(selected), len(organisms)))
    for i, prot in enumerate(selected):
        for j, org in enumerate(organisms):
            matrix[i, j] = gap_matrix.get((prot, org), 0)
    
    fig, ax = plt.subplots(figsize=(12, 14), facecolor=PALETTE['bg'])
    ax.set_facecolor(PALETTE['bg'])
    
    im = ax.imshow(matrix, cmap='RdYlGn_r', aspect='auto', vmin=0, vmax=0.5)
    
    ax.set_xticks(np.arange(len(organisms)))
    ax.set_yticks(np.arange(len(selected)))
    ax.set_xticklabels(organisms, rotation=45, ha='right', fontsize=9)
    ax.set_yticklabels(selected, fontsize=8)
    
    # Add text annotations
    for i in range(len(selected)):
        for j in range(len(organisms)):
            val = matrix[i, j]
            color = 'white' if val > 0.3 else PALETTE['text']
            ax.text(j, i, f'{val:.2f}', ha='center', va='center', fontsize=6.5, color=color)
    
    cbar = plt.colorbar(im, ax=ax, shrink=0.6, pad=0.02)
    cbar.set_label('CAI - tAI Gap', fontsize=10, color=PALETTE['text'])
    
    ax.set_title('CAI-tAI Gap: Where Codon Optimization Exceeds tRNA Capacity',
                 fontsize=13, fontweight='bold', color=PALETTE['text'], pad=16)
    ax.set_xlabel('Organism', fontsize=11, color=PALETTE['text'], labelpad=8)
    ax.set_ylabel('Protein', fontsize=11, color=PALETTE['text'], labelpad=8)
    
    plt.tight_layout()
    fig.savefig(OUTPUT_DIR / 'chart2_tai_gap_heatmap.png', dpi=300, bbox_inches='tight',
                facecolor=PALETTE['bg'])
    plt.close(fig)
    print(f"  Saved: chart2_tai_gap_heatmap.png")


def chart_3_tai_distribution_by_organism(data):
    """Violin/box plot: tAI distribution per organism."""
    tai_results = data['tai_results']
    
    org_data = {}
    for r in tai_results:
        if 'error' in r:
            continue
        org = r['organism_display']
        org_data.setdefault(org, []).append(r['tai'])
    
    organisms = sorted(org_data.keys(), key=lambda x: np.median(org_data[x]), reverse=True)
    
    fig, ax = plt.subplots(figsize=(14, 7), facecolor=PALETTE['bg'])
    ax.set_facecolor(PALETTE['bg'])
    
    positions = range(len(organisms))
    bp_data = [org_data[o] for o in organisms]
    
    bp = ax.boxplot(bp_data, positions=positions, patch_artist=True, widths=0.6,
                    medianprops=dict(color=PALETTE['text'], linewidth=2),
                    whiskerprops=dict(color=PALETTE['neutral'], linewidth=1),
                    capprops=dict(color=PALETTE['neutral'], linewidth=1),
                    flierprops=dict(marker='o', markerfacecolor=PALETTE['neutral'],
                                   markersize=4, alpha=0.5))
    
    colors = [ORGANISM_COLORS.get(o, PALETTE['primary']) for o in organisms]
    for patch, color in zip(bp['boxes'], colors):
        patch.set_facecolor(color)
        patch.set_alpha(0.6)
        patch.set_edgecolor(color)
    
    ax.set_xticklabels(organisms, rotation=30, ha='right', fontsize=10)
    ax.set_ylabel('tAI Value', fontsize=12, color=PALETTE['text'], labelpad=10)
    ax.set_title('tAI Distribution Across 59 Proteins by Organism',
                 fontsize=14, fontweight='bold', color=PALETTE['text'], pad=16)
    ax.set_ylim(0.4, 1.0)
    ax.grid(axis='y', alpha=0.15, color=PALETTE['grid'])
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.spines['left'].set_color(PALETTE['grid'])
    ax.spines['bottom'].set_color(PALETTE['grid'])
    
    plt.tight_layout()
    fig.savefig(OUTPUT_DIR / 'chart3_tai_distribution.png', dpi=300, bbox_inches='tight',
                facecolor=PALETTE['bg'])
    plt.close(fig)
    print(f"  Saved: chart3_tai_distribution.png")


def chart_4_full_pipeline_cai_tai(data):
    """Grouped bar chart: Full pipeline CAI and tAI for subset."""
    fp = data.get('full_pipeline_results', [])
    if not fp:
        print("  Skipping chart4: no full pipeline results")
        return
    
    success = [r for r in fp if r['success'] and r.get('tai', 0) > 0]
    if not success:
        print("  Skipping chart4: no successful full pipeline results with tAI")
        return
    
    # Group by organism
    org_data = {}
    for r in success:
        org = r['organism'].replace('Homo_sapiens', 'Human').replace('Escherichia_coli', 'E.coli').replace('Saccharomyces_cerevisiae', 'Yeast').replace('CHO_K1', 'CHO-K1').replace('Komagataella_phaffii', 'Pichia')
        if org not in org_data:
            org_data[org] = {'cai': [], 'tai': [], 'names': []}
        org_data[org]['cai'].append(r['cai'])
        org_data[org]['tai'].append(r['tai'])
        org_data[org]['names'].append(r['protein_name'])
    
    organisms = sorted(org_data.keys())
    
    fig, axes = plt.subplots(1, len(organisms), figsize=(5*len(organisms), 7), 
                              facecolor=PALETTE['bg'], sharey=True)
    if len(organisms) == 1:
        axes = [axes]
    
    for idx, org in enumerate(organisms):
        ax = axes[idx]
        ax.set_facecolor(PALETTE['bg'])
        
        names = org_data[org]['names']
        cais = org_data[org]['cai']
        tais = org_data[org]['tai']
        
        x = np.arange(len(names))
        width = 0.35
        
        ax.bar(x - width/2, cais, width, label='CAI',
               color=PALETTE['primary'], alpha=0.85)
        ax.bar(x + width/2, tais, width, label='tAI',
               color=PALETTE['accent'], alpha=0.85)
        
        ax.set_title(org, fontsize=12, fontweight='bold', color=PALETTE['text'])
        ax.set_xticks(x)
        ax.set_xticklabels(names, rotation=60, ha='right', fontsize=7)
        ax.set_ylim(0, 1.1)
        ax.grid(axis='y', alpha=0.15, color=PALETTE['grid'])
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)
        if idx == 0:
            ax.set_ylabel('Index Value', fontsize=11, color=PALETTE['text'])
        ax.legend(loc='upper left', fontsize=8)
    
    fig.suptitle('Full Pipeline (13-Layer) CAI and tAI Results',
                 fontsize=14, fontweight='bold', color=PALETTE['text'], y=1.02)
    plt.tight_layout()
    fig.savefig(OUTPUT_DIR / 'chart4_full_pipeline_cai_tai.png', dpi=300, bbox_inches='tight',
                facecolor=PALETTE['bg'])
    plt.close(fig)
    print(f"  Saved: chart4_full_pipeline_cai_tai.png")


def chart_5_cpg_by_organism(data):
    """Bar chart: mean CpG count by organism (CAI-optimized)."""
    cai_results = data['cai_results']
    
    # Only show eukaryotic organisms (CpG is relevant for eukaryotes)
    org_data = {}
    for r in cai_results:
        if 'error' in r or r['domain'] != 'eukaryote':
            continue
        org = r['organism_display']
        org_data.setdefault(org, []).append(r['cpg_count'])
    
    organisms = sorted(org_data.keys(), key=lambda x: np.mean(org_data[x]))
    mean_cpgs = [np.mean(org_data[o]) for o in organisms]
    std_cpgs = [np.std(org_data[o]) for o in organisms]
    
    fig, ax = plt.subplots(figsize=(14, 7), facecolor=PALETTE['bg'])
    ax.set_facecolor(PALETTE['bg'])
    
    colors = [PALETTE['warm'] if v > 30 else PALETTE['accent'] if v > 15 else PALETTE['primary']
              for v in mean_cpgs]
    
    bars = ax.barh(organisms, mean_cpgs, xerr=std_cpgs, color=colors, alpha=0.8,
                   edgecolor='white', linewidth=0.5, capsize=3)
    
    for bar, val in zip(bars, mean_cpgs):
        ax.text(bar.get_width() + 1, bar.get_y() + bar.get_height()/2,
                f'{val:.1f}', ha='left', va='center', fontsize=9, color=PALETTE['text'])
    
    ax.set_xlabel('Mean CpG Dinucleotide Count', fontsize=12, color=PALETTE['text'], labelpad=10)
    ax.set_title('CpG Dinucleotide Content by Organism (CAI-Optimized, Eukaryotes Only)',
                 fontsize=13, fontweight='bold', color=PALETTE['text'], pad=16)
    ax.grid(axis='x', alpha=0.15, color=PALETTE['grid'])
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.spines['left'].set_color(PALETTE['grid'])
    ax.spines['bottom'].set_color(PALETTE['grid'])
    
    plt.tight_layout()
    fig.savefig(OUTPUT_DIR / 'chart5_cpg_by_organism.png', dpi=300, bbox_inches='tight',
                facecolor=PALETTE['bg'])
    plt.close(fig)
    print(f"  Saved: chart5_cpg_by_organism.png")


def chart_6_gc_distribution(data):
    """Scatter: GC% vs CAI for all proteins across key organisms."""
    cai_results = data['cai_results']
    
    # Pick 5 representative organisms
    key_organisms = ['Human', 'E. coli', 'Yeast', 'CHO-K1', 'Pichia']
    
    fig, ax = plt.subplots(figsize=(12, 8), facecolor=PALETTE['bg'])
    ax.set_facecolor(PALETTE['bg'])
    
    for org in key_organisms:
        org_results = [r for r in cai_results if r['organism_display'] == org and 'error' not in r]
        if not org_results:
            continue
        gcs = [r['gc_percent'] for r in org_results]
        cais = [r['cai'] for r in org_results]
        color = ORGANISM_COLORS.get(org, PALETTE['neutral'])
        ax.scatter(gcs, cais, alpha=0.5, s=40, color=color, label=org, edgecolors='white', linewidth=0.3)
    
    ax.set_xlabel('GC Content (%)', fontsize=12, color=PALETTE['text'], labelpad=10)
    ax.set_ylabel('CAI', fontsize=12, color=PALETTE['text'], labelpad=10)
    ax.set_title('GC Content vs CAI Across Key Organisms (CAI-Optimized Sequences)',
                 fontsize=13, fontweight='bold', color=PALETTE['text'], pad=16)
    ax.legend(loc='best', fontsize=10, framealpha=0.9)
    ax.grid(alpha=0.15, color=PALETTE['grid'])
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.spines['left'].set_color(PALETTE['grid'])
    ax.spines['bottom'].set_color(PALETTE['grid'])
    
    plt.tight_layout()
    fig.savefig(OUTPUT_DIR / 'chart6_gc_vs_cai.png', dpi=300, bbox_inches='tight',
                facecolor=PALETTE['bg'])
    plt.close(fig)
    print(f"  Saved: chart6_gc_vs_cai.png")


def main():
    print("Loading benchmark data...")
    data = load_data()
    
    print("Generating charts:")
    chart_1_tai_by_organism(data)
    chart_2_tai_gap_heatmap(data)
    chart_3_tai_distribution_by_organism(data)
    chart_4_full_pipeline_cai_tai(data)
    chart_5_cpg_by_organism(data)
    chart_6_gc_distribution(data)
    
    print(f"\nAll charts saved to: {OUTPUT_DIR}")


if __name__ == '__main__':
    main()
