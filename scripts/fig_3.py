#!/usr/bin/env python3
"""
Figure 3 — Macrolide ARG abundance by AZM prophylaxis status
Publication-quality version for ERJ submission.

Layout (double-column, 175 mm wide):
  3×3 grid: 8 key macrolide genes + total macrolide RPKM (panel I)
  Each panel: box + jitter, AZM+ vs AZM−, Fisher p annotation

Data sources (real Bowtie2 mapping RPKM):
  data/arg_rpkm_mapped.tsv       — per-contig RPKM
  data/macrolide_depth_hits.tsv  — gene-to-contig assignments
  data/macrolide_rpkm_mapped.tsv — per-sample total macrolide RPKM
"""

import pandas as pd
import numpy as np
from pathlib import Path
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from scipy.stats import mannwhitneyu, fisher_exact

# ── ERJ figure specifications ────────────────────────────────────────────────
MM = 1 / 25.4
FIG_W = 175 * MM       # double column
FIG_H = 130 * MM       # 3 rows

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
BASE   = Path(__file__).parent.parent
DEPTH  = BASE / 'data' / 'macrolide_depth_hits.tsv'
OUTPDF = BASE / 'figures' / 'fig_3.pdf'
OUTPNG = BASE / 'figures' / 'fig_3.png'

# ── Colours ───────────────────────────────────────────────────────────────────
COLORS = {'AZMpos': '#C0392B', 'AZMneg': '#2980B9'}

# ── Load data ─────────────────────────────────────────────────────────────────
contigs = pd.read_csv(BASE / 'data' / 'arg_contigs.tsv', sep='\t')
depth   = pd.read_csv(DEPTH, sep='\t')
meta    = pd.read_csv(BASE / 'data' / 'sample_metadata.tsv', sep='\t')

# Join: assign real RPKM to each gene detection
gene_rpkm = depth.merge(
    contigs[['sample_id', 'contig', 'rpkm']],
    left_on=['sample_id', 'Contig id'],
    right_on=['sample_id', 'contig'],
    how='left'
)

# Total macrolide RPKM per sample (with AZM metadata)
all_samples = sorted(meta['sample_id'].tolist())
azm_map = meta.set_index('sample_id')['AZM'].to_dict()

# Per-sample per-gene RPKM (sum across contigs if gene on multiple contigs)
KEY_GENES = ['mef(A)', 'msr(D)', 'erm(B)', 'erm(X)', 'erm(F)',
             'erm(A)', 'erm(T)', 'erm(C)']   # 8 genes + total = 9 panels

gene_samp = gene_rpkm.groupby(['sample_id', 'gene'])['rpkm'].sum().reset_index()
gene_samp.columns = ['sample_id', 'gene', 'rpkm_sum']

full_idx = pd.MultiIndex.from_product([all_samples, KEY_GENES],
                                       names=['sample_id', 'gene'])
gene_full = (gene_samp
             .set_index(['sample_id', 'gene'])
             .reindex(full_idx, fill_value=0.0)
             .reset_index())
gene_full['AZM'] = gene_full['sample_id'].map(azm_map)

# ── Figure ────────────────────────────────────────────────────────────────────
fig, axes = plt.subplots(3, 3, figsize=(FIG_W, FIG_H))
axes = axes.flatten()


def draw_panel(ax, azm_pos_vals, azm_neg_vals, title, n_pos_detect, n_neg_detect,
               fisher_p, mw_p, ylabel='RPKM'):
    rng = np.random.default_rng(42)
    for jj, (vals, azm_val) in enumerate([(azm_pos_vals, 'AZMpos'),
                                           (azm_neg_vals, 'AZMneg')]):
        col = COLORS[azm_val]
        ax.boxplot(vals, positions=[jj], widths=0.45, patch_artist=True,
                   medianprops=dict(color='black', linewidth=1.2),
                   boxprops=dict(facecolor=col, alpha=0.35, linewidth=0.6),
                   whiskerprops=dict(linewidth=0.6, color='#555'),
                   capprops=dict(linewidth=0.6, color='#555'),
                   flierprops=dict(marker='', linestyle='none'))
        jitter = rng.uniform(-0.14, 0.14, len(vals))
        ax.scatter(np.full(len(vals), jj) + jitter, vals,
                   color=col, alpha=0.75, s=6, zorder=5, linewidths=0)

    p_col    = '#C0392B' if fisher_p < 0.05 else ('#E67E22' if fisher_p < 0.1 else '#888888')
    p_weight = 'bold' if fisher_p < 0.05 else 'normal'
    ax.text(0.5, 0.99,
            f'Fisher p = {fisher_p:.3f}',
            transform=ax.transAxes, ha='center', va='top',
            fontsize=5.5, color=p_col, fontweight=p_weight)

    ax.set_title(f'{title}\n(+: {n_pos_detect}/20  −: {n_neg_detect}/20)',
                 fontsize=6.5, pad=3)
    ax.set_xticks([0, 1])
    ax.set_xticklabels(['AZM+', 'AZM−'], fontsize=6)
    ax.set_ylabel(ylabel, fontsize=6)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.tick_params(axis='both', length=2, width=0.6)


# Gene panels A–H
gene_results = {}
for idx, gene in enumerate(KEY_GENES):
    ax   = axes[idx]
    gdf  = gene_full[gene_full['gene'] == gene]
    azmp = gdf[gdf['AZM'] == 'AZMpos']['rpkm_sum'].values
    azmn = gdf[gdf['AZM'] == 'AZMneg']['rpkm_sum'].values
    np_  = int((azmp > 0).sum())
    nn   = int((azmn > 0).sum())
    fe   = fisher_exact([[np_, 20 - np_], [nn, 20 - nn]])
    try:
        mw_p = mannwhitneyu(azmp, azmn, alternative='two-sided').pvalue
    except Exception:
        mw_p = np.nan
    gene_results[gene] = {'fisher_p': fe[1], 'mw_p': mw_p,
                          'n_pos': np_, 'n_neg': nn}
    draw_panel(ax, azmp, azmn, gene, np_, nn, fe[1], mw_p)

# Total macrolide panel I (index 8)
ax   = axes[8]
azmp = meta[meta['AZM'] == 'AZMpos']['total_mac_rpkm_mapped'].values
azmn = meta[meta['AZM'] == 'AZMneg']['total_mac_rpkm_mapped'].values
np_  = int((azmp > 0).sum())
nn   = int((azmn > 0).sum())
fe   = fisher_exact([[np_, 20 - np_], [nn, 20 - nn]])
mw   = mannwhitneyu(azmp, azmn, alternative='two-sided')

ax_rng = np.random.default_rng(99)
for jj, (vals, azm_val) in enumerate([(azmp, 'AZMpos'), (azmn, 'AZMneg')]):
    col = COLORS[azm_val]
    ax.boxplot(vals, positions=[jj], widths=0.45, patch_artist=True,
               medianprops=dict(color='black', linewidth=1.2),
               boxprops=dict(facecolor=col, alpha=0.35, linewidth=0.6),
               whiskerprops=dict(linewidth=0.6, color='#555'),
               capprops=dict(linewidth=0.6, color='#555'),
               flierprops=dict(marker='', linestyle='none'))
    jitter = ax_rng.uniform(-0.14, 0.14, len(vals))
    ax.scatter(np.full(len(vals), jj) + jitter, vals,
               color=col, alpha=0.75, s=6, zorder=5, linewidths=0)

ax.text(0.5, 0.99, f'MW p = {mw.pvalue:.3f}',
        transform=ax.transAxes, ha='center', va='top',
        fontsize=5.5, color='#888888')
ax.set_title(f'Total macrolide\n(all genes)', fontsize=6.5, pad=3)
ax.set_xticks([0, 1])
ax.set_xticklabels(['AZM+', 'AZM−'], fontsize=6)
ax.set_ylabel('Total RPKM', fontsize=6)
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
ax.tick_params(axis='both', length=2, width=0.6)

# ─── Panel labels ─────────────────────────────────────────────────────────────
for i, ax in enumerate(axes):
    label = chr(65 + i)   # A, B, C, …
    ax.text(-0.14, 1.06, label, transform=ax.transAxes,
            fontsize=9, fontweight='bold', va='top')

fig.suptitle(
    'Macrolide ARG abundance by azithromycin prophylaxis status\n'
    '(assembly-based, RGI + AMRFinderPlus, n = 40 CF patients)',
    fontsize=7, y=0.995
)
plt.tight_layout(rect=[0, 0, 1, 0.975])

# ─── Save ─────────────────────────────────────────────────────────────────────
fig.savefig(OUTPDF, bbox_inches='tight', dpi=300)
fig.savefig(OUTPNG, bbox_inches='tight', dpi=300)
plt.close(fig)
print(f'Saved:\n  {OUTPDF}\n  {OUTPNG}')
print()
print('=== Gene-level results (real mapping RPKM) ===')
for gene, res in gene_results.items():
    print(f"  {gene:12s}  AZM+ {res['n_pos']:2d}/20  AZM- {res['n_neg']:2d}/20"
          f"  Fisher p={res['fisher_p']:.3f}  MW p={res['mw_p']:.3f}")
print(f"  {'Total':12s}  MW p={mw.pvalue:.3f}")
