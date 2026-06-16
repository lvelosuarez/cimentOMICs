#!/usr/bin/env python3
"""
Figure 1 — PA pulmotype drives bronchopulmonary ARG richness
Publication-quality version for ERJ submission.

Layout (double-column, 175 mm wide):
  Panel A: Stacked bar — ARG richness by organism group, samples sorted by PA%
  Panel B: Scatter — PA% vs total ARG richness (Spearman correlation)
  Panel C: Scatter — PA% vs non-PA ARG richness (Spearman correlation)
"""

import pandas as pd
import numpy as np
from pathlib import Path
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.lines import Line2D
from scipy.stats import spearmanr

# ── ERJ figure specifications ────────────────────────────────────────────────
MM = 1 / 25.4          # mm → inches
FIG_W = 175 * MM       # double column
FIG_H = 75 * MM        # height (extra space for organism legend below)

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
    'pdf.fonttype':     42,   # embed fonts as TrueType
    'ps.fonttype':      42,
})

# ── Paths ────────────────────────────────────────────────────────────────────
BASE  = Path(__file__).parent.parent
DATA  = BASE / 'data' / 'assembly_arg_richness_taxonomy.tsv'
OUTPDF = BASE / 'figures' / 'fig_2.pdf'
OUTPNG = BASE / 'figures' / 'fig_2.png'

# ── Colour palette (Okabe-Ito, colorblind-safe) ──────────────────────────────
# Organism groups — consistent order throughout
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

ORG_COLORS = {
    'Pseudomonas aeruginosa':           '#0072B2',  # deep blue
    'Achromobacter spp.':               '#E69F00',  # amber
    'Stenotrophomonas spp.':            '#009E73',  # teal
    'Haemophilus spp.':                 '#56B4E9',  # sky blue
    'Streptococcus / oral Firmicutes':  '#D55E00',  # vermillion
    'Oral anaerobes':                   '#CC79A7',  # reddish purple
    'Rothia / Actinomyces':             '#F0E442',  # yellow
    'Neisseria spp.':                   '#7570B3',  # purple
    'Staphylococcus aureus':            '#1B9E77',  # green
    'Burkholderiales (unclassified)':   '#A6761D',  # brown
    'Enterobacteria / other':           '#E7298A',  # pink
    'Other/Unclassified':               '#BDBDBD',  # light grey
}

AZM_COLORS = {'AZMpos': '#C0392B', 'AZMneg': '#2980B9'}  # red / blue

# ── Load data ─────────────────────────────────────────────────────────────────
df = pd.read_csv(DATA, sep='\t')
df = df.sort_values('PA_pct').reset_index(drop=True)

# ── Figure layout ─────────────────────────────────────────────────────────────
fig = plt.figure(figsize=(FIG_W, FIG_H))
# Panels: A wider (stacked bar), B and C equal scatter
gs = fig.add_gridspec(
    1, 3,
    width_ratios=[2.2, 1, 1],
    left=0.07, right=0.98,
    top=0.92, bottom=0.22,
    wspace=0.38,
)
axA = fig.add_subplot(gs[0])
axB = fig.add_subplot(gs[1])
axC = fig.add_subplot(gs[2])

# ─── Panel A: stacked bar ─────────────────────────────────────────────────────
x = np.arange(len(df))
bottom = np.zeros(len(df))
for org in ORG_ORDER:
    if org not in df.columns:
        continue
    vals = df[org].values.astype(float)
    axA.bar(x, vals, bottom=bottom,
            color=ORG_COLORS[org], linewidth=0, width=0.9)
    bottom += vals

# PA% overlay line (secondary axis)
ax_twin = axA.twinx()
ax_twin.plot(x, df['PA_pct'].values, color='black', linewidth=0.9,
             linestyle='-', zorder=5)
ax_twin.set_ylabel('P. aeruginosa\nrelative abundance (%)', fontsize=6,
                   labelpad=3)
ax_twin.tick_params(axis='y', labelsize=6, length=2.5, width=0.6)
ax_twin.set_ylim(0, df['PA_pct'].max() * 1.15)
ax_twin.spines['top'].set_visible(False)

# AZM markers below x-axis
for i, row in df.iterrows():
    col = AZM_COLORS[row['AZM']]
    axA.plot(i, -1.5, 's', color=col, markersize=2.5, clip_on=False)

axA.set_xlim(-0.7, len(df) - 0.3)
axA.set_ylim(0, df['n_total'].max() * 1.05)
axA.set_xticks([])
axA.set_xlabel('Patients (sorted by P. aeruginosa abundance)', fontsize=7,
               labelpad=8)
axA.set_ylabel('ARG richness (no. unique ARO accessions)', fontsize=7)
axA.spines['top'].set_visible(False)
axA.spines['right'].set_visible(False)

# Organism legend — placed below panel A in the figure coordinate space
handles_org = [mpatches.Patch(color=ORG_COLORS[o], label=o)
               for o in ORG_ORDER if o in df.columns]
leg_org = fig.legend(
    handles=handles_org,
    loc='lower left',
    bbox_to_anchor=(0.065, -0.01),
    bbox_transform=fig.transFigure,
    ncol=3, fontsize=5, frameon=False,
    handlelength=1, handleheight=0.8,
    columnspacing=0.6, labelspacing=0.25, borderpad=0,
)

# AZM legend — small coloured squares below x-axis tick area
azm_handles = [
    Line2D([0], [0], marker='s', color='w', markerfacecolor=AZM_COLORS['AZMpos'],
           markersize=4, label='AZM+'),
    Line2D([0], [0], marker='s', color='w', markerfacecolor=AZM_COLORS['AZMneg'],
           markersize=4, label='AZM−'),
]
axA.legend(handles=azm_handles, loc='upper left', fontsize=6,
           frameon=False, ncol=2, handlelength=0.8,
           borderpad=0, labelspacing=0.2)

# ─── Panels B & C: scatter plots ─────────────────────────────────────────────
def scatter_panel(ax, x_vals, y_vals, azm_col, xlabel, ylabel, rho, pval):
    for azm_val, color in AZM_COLORS.items():
        mask = azm_col == azm_val
        ax.scatter(x_vals[mask], y_vals[mask],
                   c=color, s=12, alpha=0.8, linewidths=0,
                   label='AZM+' if azm_val == 'AZMpos' else 'AZM−', zorder=3)
    # Regression line (all samples)
    m, b = np.polyfit(x_vals, y_vals, 1)
    xr = np.linspace(x_vals.min(), x_vals.max(), 100)
    ax.plot(xr, m * xr + b, color='#444444', linewidth=0.8, linestyle='--', zorder=2)
    # Stats annotation
    p_str = 'p < 0.001' if pval < 0.001 else f'p = {pval:.3f}'
    ax.text(0.97, 0.05, f'ρ = {rho:.2f}\n{p_str}',
            transform=ax.transAxes, ha='right', va='bottom', fontsize=6,
            fontstyle='italic')
    ax.set_xlabel(xlabel, fontsize=7)
    ax.set_ylabel(ylabel, fontsize=7)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.tick_params(axis='both', length=2.5, width=0.6)

pa_pct  = df['PA_pct'].values
n_total = df['n_total'].values
n_nonPA = df['n_nonPA'].values
azm_col = df['AZM']

rho_total, p_total = spearmanr(pa_pct, n_total)
rho_nonPA, p_nonPA = spearmanr(pa_pct, n_nonPA)

scatter_panel(axB, pa_pct, n_total, azm_col,
              'P. aeruginosa\nabundance (%)',
              'Total ARG richness', rho_total, p_total)

scatter_panel(axC, pa_pct, n_nonPA, azm_col,
              'P. aeruginosa\nabundance (%)',
              'Non-PA ARG richness', rho_nonPA, p_nonPA)

# Shared legend for B/C
handles_sc = [
    Line2D([0], [0], marker='o', color='w', markerfacecolor=AZM_COLORS['AZMpos'],
           markersize=4, label='AZM+'),
    Line2D([0], [0], marker='o', color='w', markerfacecolor=AZM_COLORS['AZMneg'],
           markersize=4, label='AZM−'),
]
axC.legend(handles=handles_sc, fontsize=6, frameon=False,
           loc='upper left', handlelength=0.8)

# ─── Panel labels ─────────────────────────────────────────────────────────────
for ax, label in [(axA, 'A'), (axB, 'B'), (axC, 'C')]:
    ax.text(-0.10, 1.06, label, transform=ax.transAxes,
            fontsize=9, fontweight='bold', va='top')

# ─── Save ─────────────────────────────────────────────────────────────────────
fig.savefig(OUTPDF, bbox_inches='tight', dpi=300)
fig.savefig(OUTPNG, bbox_inches='tight', dpi=300)
plt.close(fig)
print(f'Saved:\n  {OUTPDF}\n  {OUTPNG}')
