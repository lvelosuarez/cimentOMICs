#!/usr/bin/env python3
"""
Figure 4 — Mobile genetic element context of CF bronchopulmonary ARGs
Publication-quality version for ERJ submission.

Layout (double-column, 175 mm wide):
  Panel A: Stacked bar — mobile vs chromosomal ARG contigs per organism group
  Panel B: Box+jitter — mobile ARG count per sample by AZM group (p=0.913)
  Panel C: Pie — overall mobile (16.3%) vs chromosomal (83.7%)
  Panel D: Pie — macrolide ARG contigs (87.9% chromosomal)
  Panel E: Bar — mobile macrolide ARG gene breakdown

Data sources:
  data/arg_mge_colocalisation.tsv    — all ARG contigs with on_plasmid flag
  data/macrolide_depth_hits.tsv      — gene-to-contig for macrolide ARGs
  data/genomad_mobile_per_sample.tsv — per-sample mobile ARG counts
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
FIG_H = 110 * MM

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
OUTPDF  = BASE / 'figures' / 'fig_4.pdf'
OUTPNG  = BASE / 'figures' / 'fig_4.png'

# ── Colours ───────────────────────────────────────────────────────────────────
COL_MOB   = '#E69F00'   # amber — mobile
COL_CHROM = '#0072B2'   # blue — chromosomal
AZM_COLORS = {'AZMpos': '#C0392B', 'AZMneg': '#2980B9'}

ORG_ORDER = [
    'Pseudomonas aeruginosa',
    'Achromobacter spp.',
    'Stenotrophomonas spp.',
    'Haemophilus spp.',
    'Streptococcus / oral Firmicutes',
    'Oral anaerobes',
    'Rothia / Actinomyces',
    'Neisseria spp.',
    'Staphylococcus aureus',
    'Burkholderiales (unclassified)',
    'Enterobacteria / other',
    'Other/Unclassified',
]

# ── Load data ─────────────────────────────────────────────────────────────────
# All ARG contigs with MGE classification
args = pd.read_csv(BASE / 'data' / 'arg_contigs.tsv', sep='\t')

# Per-sample mobile ARG counts with AZM metadata
mobile = pd.read_csv(BASE / 'data' / 'sample_metadata.tsv', sep='\t')

# Macrolide contig list — join with args to get on_plasmid flag
depth = pd.read_csv(BASE / 'data' / 'macrolide_depth_hits.tsv', sep='\t')
mac_mge = depth.merge(
    args[['sample_id', 'contig', 'on_plasmid']],
    left_on=['sample_id', 'Contig id'],
    right_on=['sample_id', 'contig'],
    how='left'
)

# ── Figure ────────────────────────────────────────────────────────────────────
fig = plt.figure(figsize=(FIG_W, FIG_H))
gs = fig.add_gridspec(
    2, 3,
    left=0.08, right=0.97,
    top=0.92, bottom=0.12,
    hspace=0.50, wspace=0.42,
    width_ratios=[1.6, 1, 1],
)
axA = fig.add_subplot(gs[:, 0])   # tall left panel
axB = fig.add_subplot(gs[0, 1])
axC = fig.add_subplot(gs[0, 2])
axD = fig.add_subplot(gs[1, 1])
axE = fig.add_subplot(gs[1, 2])

# ─── Panel A: stacked bar per organism ───────────────────────────────────────
org_mob   = args[args['on_plasmid']]['broad_group'].value_counts()
org_chrom = args[~args['on_plasmid']]['broad_group'].value_counts()

org_total = org_mob.add(org_chrom, fill_value=0)
orgs = [o for o in ORG_ORDER if o in org_total.index][::-1]

y = np.arange(len(orgs))
chrom_vals = [org_chrom.get(o, 0) for o in orgs]
mob_vals   = [org_mob.get(o, 0) for o in orgs]

axA.barh(y, chrom_vals, color=COL_CHROM, height=0.65, label='Chromosomal')
axA.barh(y, mob_vals,   left=chrom_vals, color=COL_MOB, height=0.65, label='Plasmid-associated')

for i, (c, m) in enumerate(zip(chrom_vals, mob_vals)):
    total = c + m
    if total > 0 and m > 0:
        pct = 100 * m / total
        axA.text(total + 1, i, f'{pct:.0f}%', va='center', fontsize=5, color='#444')

axA.set_yticks(y)
axA.set_yticklabels([o.replace(' / ', '/') for o in orgs], fontsize=6)
axA.set_xlabel('No. ARG-carrying contigs', fontsize=7)
axA.spines['top'].set_visible(False)
axA.spines['right'].set_visible(False)
axA.tick_params(axis='both', length=2.5, width=0.6)
axA.legend(fontsize=6, frameon=False, loc='lower right',
           handlelength=1, handleheight=0.8)

# ─── Panel B: AZM comparison — mobile ARG count ──────────────────────────────
rng = np.random.default_rng(42)
for jj, (azm_val, color) in enumerate(AZM_COLORS.items()):
    vals = mobile[mobile['AZM'] == azm_val]['n_mobile_arg'].values
    axB.boxplot(vals, positions=[jj], widths=0.45, patch_artist=True,
                medianprops=dict(color='black', linewidth=1.2),
                boxprops=dict(facecolor=color, alpha=0.35, linewidth=0.6),
                whiskerprops=dict(linewidth=0.6, color='#555'),
                capprops=dict(linewidth=0.6, color='#555'),
                flierprops=dict(marker='', linestyle='none'))
    jitter = rng.uniform(-0.14, 0.14, len(vals))
    axB.scatter(np.full(len(vals), jj) + jitter, vals,
                color=color, alpha=0.75, s=8, zorder=5, linewidths=0)

azmp = mobile[mobile['AZM'] == 'AZMpos']['n_mobile_arg']
azmn = mobile[mobile['AZM'] == 'AZMneg']['n_mobile_arg']
mw   = mannwhitneyu(azmp, azmn, alternative='two-sided')
axB.text(0.5, 0.98, f'MW p = {mw.pvalue:.3f}',
         transform=axB.transAxes, ha='center', va='top',
         fontsize=5.5, color='#888888')
axB.set_xticks([0, 1])
axB.set_xticklabels(['AZM+', 'AZM−'], fontsize=6)
axB.set_ylabel('Mobile ARG contigs per patient', fontsize=7)
axB.spines['top'].set_visible(False)
axB.spines['right'].set_visible(False)
axB.tick_params(axis='both', length=2.5, width=0.6)

# ─── Panel C: Overall pie — mobile vs chromosomal ────────────────────────────
n_mob   = int(args['on_plasmid'].sum())
n_chrom = int((~args['on_plasmid']).sum())
wedge_colors = [COL_CHROM, COL_MOB]
wedges, texts, autotexts = axC.pie(
    [n_chrom, n_mob],
    labels=[f'Chromosomal\n({n_chrom})', f'Plasmid-\nassociated\n({n_mob})'],
    colors=wedge_colors,
    autopct='%1.1f%%',
    startangle=90,
    textprops={'fontsize': 6},
    wedgeprops={'linewidth': 0.5, 'edgecolor': 'white'},
    pctdistance=0.65,
)
for at in autotexts:
    at.set_fontsize(6)
    at.set_fontweight('bold')
    at.set_color('white')
axC.set_title('All ARG-carrying\ncontigs', fontsize=7, pad=3)

# ─── Panel D: Macrolide ARG contigs — mobile vs chromosomal ──────────────────
n_mac_mob   = int(mac_mge['on_plasmid'].sum())
n_mac_chrom = int((~mac_mge['on_plasmid']).sum())
wedges2, texts2, autotexts2 = axD.pie(
    [n_mac_chrom, n_mac_mob],
    labels=[f'Chromosomal\n({n_mac_chrom})', f'Plasmid-\nassociated\n({n_mac_mob})'],
    colors=wedge_colors,
    autopct='%1.1f%%',
    startangle=90,
    textprops={'fontsize': 6},
    wedgeprops={'linewidth': 0.5, 'edgecolor': 'white'},
    pctdistance=0.65,
)
for at in autotexts2:
    at.set_fontsize(6)
    at.set_fontweight('bold')
    at.set_color('white')
axD.set_title('Macrolide ARG\ncontigs', fontsize=7, pad=3)

# ─── Panel E: Mobile macrolide gene breakdown ─────────────────────────────────
mob_mac = mac_mge[mac_mge['on_plasmid'] == True].copy()
gene_counts = mob_mac['gene'].value_counts()
axE.barh(np.arange(len(gene_counts)), gene_counts.values,
         color=COL_MOB, height=0.65)
axE.set_yticks(np.arange(len(gene_counts)))
axE.set_yticklabels(gene_counts.index, fontsize=6, fontstyle='italic')
axE.set_xlabel('No. plasmid-associated\nmacrolide ARG contigs', fontsize=6)
axE.set_title('Mobile macrolide ARGs', fontsize=7, pad=3)
axE.spines['top'].set_visible(False)
axE.spines['right'].set_visible(False)
axE.tick_params(axis='both', length=2.5, width=0.6)
axE.xaxis.set_major_locator(plt.MaxNLocator(integer=True))

# ─── Panel labels ─────────────────────────────────────────────────────────────
for ax, label, xoff in [(axA, 'A', -0.06), (axB, 'B', -0.18),
                         (axC, 'C', -0.12), (axD, 'D', -0.18), (axE, 'E', -0.18)]:
    ax.text(xoff, 1.06, label, transform=ax.transAxes,
            fontsize=9, fontweight='bold', va='top')

# ─── Save ─────────────────────────────────────────────────────────────────────
fig.savefig(OUTPDF, bbox_inches='tight', dpi=300)
fig.savefig(OUTPNG, bbox_inches='tight', dpi=300)
plt.close(fig)
print(f'Saved:\n  {OUTPDF}\n  {OUTPNG}')
print(f'\nKey numbers:')
print(f'  Overall: {n_mob}/{n_mob+n_chrom} mobile ({100*n_mob/(n_mob+n_chrom):.1f}%)')
print(f'  Macrolide: {n_mac_mob}/{n_mac_mob+n_mac_chrom} mobile '
      f'({100*n_mac_mob/(n_mac_mob+n_mac_chrom):.1f}%)')
print(f'\nMobile macrolide gene breakdown:')
print(gene_counts.to_string())
