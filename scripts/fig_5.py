#!/usr/bin/env python3
"""
Figure 5 — CFTR modulator use is not associated with mobile ARG enrichment
Publication-quality version for ERJ submission.

Layout (double-column, 175 mm wide):
  Panel A: Box+jitter — total ARG richness CFTR+ vs CFTR−
  Panel B: Box+jitter — PA% (sylph) CFTR+ vs CFTR−
  Panel C: Box+jitter — total macrolide RPKM CFTR+ vs CFTR−
  Panel D: Box+jitter — mobile ARG contig count CFTR+ vs CFTR−

Data sources:
  data/cftr_status.tsv                    — CFTR modulator status per sample
  data/assembly_arg_richness_taxonomy.tsv  — total ARG richness per sample
  data/sylph_ciment_summary.tsv            — PA_sylph_pct per sample
  data/macrolide_rpkm_mapped.tsv           — total macrolide RPKM (mapped)
  data/genomad_mobile_per_sample.tsv       — mobile ARG contig count

CFTR status derived from:
  data/cftr_status.tsv (pre-computed; source: clinical data not shared per ethics)
  n = 21 CFTR+ (ivacaftor/orkambi); n = 19 CFTR−
"""

import pandas as pd
import numpy as np
from pathlib import Path
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from scipy.stats import mannwhitneyu

# ── ERJ figure specifications ─────────────────────────────────────────────────
MM = 1 / 25.4
FIG_W = 175 * MM
FIG_H = 80 * MM    # 1 row of 4 panels

matplotlib.rcParams.update({
    'font.family':      'sans-serif',
    'font.sans-serif':  ['Arial', 'Helvetica', 'DejaVu Sans'],
    'font.size':        7,
    'axes.labelsize':   7,
    'axes.titlesize':   7,
    'xtick.labelsize':  6,
    'ytick.labelsize':  6,
    'legend.fontsize':  6,
    'axes.linewidth':   0.6,
    'xtick.major.width': 0.6,
    'ytick.major.width': 0.6,
    'xtick.major.size':  2.5,
    'ytick.major.size':  2.5,
    'pdf.fonttype':     42,
    'ps.fonttype':      42,
})

# ── Paths ─────────────────────────────────────────────────────────────────────
BASE    = Path(__file__).parent.parent
OUTPDF  = BASE / 'figures' / 'fig_5.pdf'
OUTPNG  = BASE / 'figures' / 'fig_5.png'

# ── Colours ───────────────────────────────────────────────────────────────────
COL_POS = '#27AE60'   # green — CFTR+
COL_NEG = '#E91E8C'   # pink  — CFTR−

# ── Load data ─────────────────────────────────────────────────────────────────
meta = pd.read_csv(BASE / 'data' / 'sample_metadata.tsv', sep='\t')
# ARG richness: unique ARO accessions per sample (RGI + AMRFinderPlus)
aro_raw  = pd.read_csv(BASE / 'data' / 'merged_aro_per_sample.tsv', sep='\t')
richness = aro_raw.groupby('sample_id')['ARO_id'].nunique().reset_index()
richness.columns = ['sample_id', 'n_aro']

# Build CFTR metadata
cftr_map = dict(zip(meta['sample_id'], meta['cftr_pos'].map({True: 'CFTR+', False: 'CFTR−'})))
n_pos = sum(v == 'CFTR+' for v in cftr_map.values())
n_neg = sum(v == 'CFTR−' for v in cftr_map.values())
print(f'CFTR+: n={n_pos}, CFTR−: n={n_neg}')

richness['CFTR'] = richness['sample_id'].map(cftr_map)
meta['CFTR']     = meta['sample_id'].map(cftr_map)

# ── Helper: draw box+jitter panel ─────────────────────────────────────────────
def draw_panel(ax, pos_vals, neg_vals, ylabel, seed=42):
    rng = np.random.default_rng(seed)
    for jj, (vals, color) in enumerate([(pos_vals, COL_POS), (neg_vals, COL_NEG)]):
        ax.boxplot(vals, positions=[jj], widths=0.45, patch_artist=True,
                   medianprops=dict(color='black', linewidth=1.2),
                   boxprops=dict(facecolor=color, alpha=0.35, linewidth=0.6),
                   whiskerprops=dict(linewidth=0.6, color='#555'),
                   capprops=dict(linewidth=0.6, color='#555'),
                   flierprops=dict(marker='', linestyle='none'))
        jitter = rng.uniform(-0.14, 0.14, len(vals))
        ax.scatter(np.full(len(vals), jj) + jitter, vals,
                   color=color, alpha=0.75, s=6, zorder=5, linewidths=0)

    mw = mannwhitneyu(pos_vals, neg_vals, alternative='two-sided')
    p_col = '#C0392B' if mw.pvalue < 0.05 else ('#E67E22' if mw.pvalue < 0.10 else '#888888')
    ax.text(0.5, 0.99, f'MW p = {mw.pvalue:.3f}',
            transform=ax.transAxes, ha='center', va='top',
            fontsize=5.5, color=p_col)
    ax.set_xticks([0, 1])
    ax.set_xticklabels([f'CFTR+\n(n={n_pos})', f'CFTR−\n(n={n_neg})'], fontsize=6)
    ax.set_ylabel(ylabel, fontsize=6)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.tick_params(axis='both', length=2, width=0.6)
    return mw.pvalue

# ── Figure ────────────────────────────────────────────────────────────────────
fig, axes = plt.subplots(1, 4, figsize=(FIG_W, FIG_H))

# Panel A: total ARG richness (unique ARO accessions)
pos_a = richness[richness['CFTR'] == 'CFTR+']['n_aro'].values
neg_a = richness[richness['CFTR'] == 'CFTR−']['n_aro'].values
p_a = draw_panel(axes[0], pos_a, neg_a, 'Total ARG richness\n(unique AROs)')
axes[0].set_title('ARG richness', fontsize=7, pad=3)

# Panel B: PA abundance
pos_b = meta[meta['CFTR'] == 'CFTR+']['PA_sylph_pct'].values
neg_b = meta[meta['CFTR'] == 'CFTR−']['PA_sylph_pct'].values
p_b = draw_panel(axes[1], pos_b, neg_b, 'PA abundance (%)')
axes[1].set_title('P. aeruginosa\nabundance (PA%)', fontsize=7, pad=3)

# Panel C: total macrolide RPKM
pos_c = meta[meta['CFTR'] == 'CFTR+']['total_mac_rpkm_mapped'].values
neg_c = meta[meta['CFTR'] == 'CFTR−']['total_mac_rpkm_mapped'].values
p_c = draw_panel(axes[2], pos_c, neg_c, 'Total macrolide\nRPKM')
axes[2].set_title('Macrolide ARG\nabundance', fontsize=7, pad=3)

# Panel D: mobile ARG contigs
pos_d = meta[meta['CFTR'] == 'CFTR+']['n_mobile_arg'].values
neg_d = meta[meta['CFTR'] == 'CFTR−']['n_mobile_arg'].values
p_d = draw_panel(axes[3], pos_d, neg_d, 'Mobile ARG\ncontigs per patient')
axes[3].set_title('Mobile ARG\ncontigs', fontsize=7, pad=3)

# Panel labels
for i, ax in enumerate(axes):
    ax.text(-0.18, 1.06, chr(65 + i), transform=ax.transAxes,
            fontsize=9, fontweight='bold', va='top')

plt.tight_layout()

# ─── Save ─────────────────────────────────────────────────────────────────────
fig.savefig(OUTPDF, bbox_inches='tight', dpi=300)
fig.savefig(OUTPNG, bbox_inches='tight', dpi=300)
plt.close(fig)
print(f'Saved:\n  {OUTPDF}\n  {OUTPNG}')
print(f'\nKey results:')
print(f'  Panel A (ARG richness):  CFTR+ {np.median(pos_a):.1f} vs CFTR− {np.median(neg_a):.1f}, p={p_a:.3f}')
print(f'  Panel B (PA%):           CFTR+ {np.median(pos_b):.1f} vs CFTR− {np.median(neg_b):.1f}, p={p_b:.3f}')
print(f'  Panel C (macrolide):     CFTR+ {np.median(pos_c):.1f} vs CFTR− {np.median(neg_c):.1f}, p={p_c:.3f}')
print(f'  Panel D (mobile ARGs):   CFTR+ {np.median(pos_d):.1f} vs CFTR− {np.median(neg_d):.1f}, p={p_d:.3f}')
