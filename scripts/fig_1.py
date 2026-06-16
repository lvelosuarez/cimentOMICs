#!/usr/bin/env python3
"""
Figure — Full bronchopulmonary resistome overview
Publication-quality for ERJ submission.

Layout (double-column, 175 mm wide):
  Row 0 (thin): PA% strip above heatmap
  Row 1 (heatmap): drug class × patient, sorted by PA%
  Row 2: lollipop top AROs | stacked group bar | mechanism donut
"""

import pandas as pd
import numpy as np
from pathlib import Path
import re
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import matplotlib.gridspec as gridspec
from matplotlib.colors import LinearSegmentedColormap, PowerNorm
from matplotlib.lines import Line2D

# ── ERJ specs ──────────────────────────────────────────────────────────────
MM = 1 / 25.4
FIG_W = 175 * MM
FIG_H = 155 * MM

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

# ── Paths ──────────────────────────────────────────────────────────────────
BASE   = Path(__file__).parent.parent
OUTPDF = BASE / 'figures' / 'fig_1.pdf'
OUTPNG = BASE / 'figures' / 'fig_1.png'

# ── Colours ────────────────────────────────────────────────────────────────
CLASS_COLORS = {
    'beta-lactam':      '#0072B2',
    'macrolide/MLSB':   '#D55E00',
    'fluoroquinolone':  '#009E73',
    'tetracycline':     '#E69F00',
    'aminoglycoside':   '#56B4E9',
    'phenicol':         '#CC79A7',
    'triclosan':        '#F0E442',
    'fosfomycin':       '#7570B3',
    'peptide':          '#1B9E77',
    'sulfonamide':      '#A6761D',
    'other':            '#BDBDBD',
}
CLASS_LABELS = {
    'beta-lactam': 'Beta-lactam', 'macrolide/MLSB': 'Macrolide/MLSB',
    'fluoroquinolone': 'Fluoroquinolone', 'tetracycline': 'Tetracycline',
    'aminoglycoside': 'Aminoglycoside', 'phenicol': 'Phenicol',
    'triclosan': 'Triclosan', 'fosfomycin': 'Fosfomycin',
    'peptide': 'Peptide', 'sulfonamide': 'Sulfonamide', 'other': 'Other',
}
CLASS_ORDER = ['beta-lactam','macrolide/MLSB','fluoroquinolone','tetracycline',
               'aminoglycoside','fosfomycin','peptide','phenicol','triclosan',
               'sulfonamide','other']
AZM_COLORS = {'AZMpos': '#C0392B', 'AZMneg': '#2980B9'}
MECH_COLORS = {
    'Efflux':              '#0072B2',
    'Inactivation':        '#D55E00',
    'Target alteration':   '#009E73',
    'Target alt.+efflux':  '#56B4E9',
    'Target protection':   '#E69F00',
    'Target replacement':  '#CC79A7',
    'Other':               '#BDBDBD',
}

# ── Load data ──────────────────────────────────────────────────────────────
BROAD_MAP = {
    'beta-lactam':     ['penam','cephalosporin','carbapenem','monobactam','cephamycin','penem','BETA-LACTAM','beta-lactam'],
    'macrolide/MLSB':  ['macrolide','lincosamide','streptogramin','MACROLIDE'],
    'aminoglycoside':  ['aminoglycoside','AMINOGLYCOSIDE'],
    'fluoroquinolone': ['fluoroquinolone','FLUOROQUINOLONE'],
    'tetracycline':    ['tetracycline','TETRACYCLINE'],
    'phenicol':        ['phenicol','PHENICOL','chloramphenicol'],
    'trimethoprim':    ['diaminopyrimidine','trimethoprim'],
    'fosfomycin':      ['fosfomycin'],
    'peptide':         ['peptide antibiotic'],
    'triclosan':       ['triclosan'],
    'sulfonamide':     ['sulfonamide'],
}
def assign_broad(s):
    if pd.isna(s): return 'other'
    sl = s.lower()
    for cat, kws in BROAD_MAP.items():
        for kw in kws:
            if kw.lower() in sl: return cat
    return 'other'

df    = pd.read_csv(BASE / 'data' / 'merged_aro_per_sample.tsv', sep='\t')
pivot = pd.read_csv(BASE / 'data' / 'resistome_class_per_sample.tsv', sep='\t').set_index('sample_id')
prev  = pd.read_csv(BASE / 'data' / 'aro_prevalence.tsv', sep='\t')
meta  = pd.read_csv(BASE / 'data' / 'sample_metadata.tsv', sep='\t').set_index('sample_id')
df['broad_class'] = df['drug_class'].apply(assign_broad)

co_samples = sorted(
    [s for s in pivot.index if s.startswith('CO')],
    key=lambda s: meta.loc[s, 'PA_sylph_pct'] if s in meta.index else 0
)
pivot_co = pivot.loc[co_samples]
classes_in_data = [c for c in CLASS_ORDER if c in pivot_co.columns]

# ── Figure layout ──────────────────────────────────────────────────────────
fig = plt.figure(figsize=(FIG_W, FIG_H))
gs = gridspec.GridSpec(
    3, 3,
    figure=fig,
    left=0.10, right=0.97,
    top=0.94, bottom=0.07,
    hspace=0.55, wspace=0.45,
    height_ratios=[0.18, 0.82, 1.0],
    width_ratios=[1.5, 1.1, 0.9],
)
ax_pa   = fig.add_subplot(gs[0, 0])          # PA% strip
ax_heat = fig.add_subplot(gs[1, 0])          # heatmap
ax_lol  = fig.add_subplot(gs[:2, 1])         # lollipop (spans rows 0-1)
ax_mech = fig.add_subplot(gs[:2, 2])         # mechanism donut (spans rows 0-1)
ax_bar  = fig.add_subplot(gs[2, :])          # group stacked bar (full width)

# ─── PA% strip ────────────────────────────────────────────────────────────
pa_vals = [meta.loc[s, 'PA_sylph_pct'] if s in meta.index else 0 for s in co_samples]
x = np.arange(len(co_samples))
ax_pa.fill_between(x, pa_vals, color='#0072B2', alpha=0.35, linewidth=0)
ax_pa.plot(x, pa_vals, color='#0072B2', linewidth=0.9)
ax_pa.set_xlim(-0.5, len(co_samples) - 0.5)
ax_pa.set_ylim(0, max(pa_vals) * 1.15)
ax_pa.set_xticks([])
ax_pa.set_yticks([0, 50, 100])
ax_pa.tick_params(axis='y', labelsize=5, length=2, width=0.5, pad=1)
ax_pa.set_ylabel('PA%', fontsize=5.5, labelpad=2)
ax_pa.spines['top'].set_visible(False)
ax_pa.spines['right'].set_visible(False)
ax_pa.spines['bottom'].set_visible(False)

# AZM strip inside PA panel
for i, s in enumerate(co_samples):
    azm = meta.loc[s, 'AZM'] if s in meta.index else 'AZMneg'
    ax_pa.plot(i, -8, 's', color=AZM_COLORS.get(azm, '#aaa'),
               markersize=2.5, clip_on=False, transform=ax_pa.transData)

azm_h = [Line2D([0],[0],marker='s',color='w',markerfacecolor=AZM_COLORS['AZMpos'],markersize=4,label='AZM+'),
         Line2D([0],[0],marker='s',color='w',markerfacecolor=AZM_COLORS['AZMneg'],markersize=4,label='AZM−')]
ax_pa.legend(handles=azm_h, fontsize=5.5, frameon=False, loc='upper left',
             handlelength=0.8, ncol=2, borderpad=0)

# ─── Heatmap ───────────────────────────────────────────────────────────────
heat_data = pivot_co[classes_in_data].T
cmap = LinearSegmentedColormap.from_list('res', ['#FFFFFF', '#0072B2'], N=256)
norm = PowerNorm(gamma=0.4, vmin=0, vmax=heat_data.values.max())
im = ax_heat.imshow(heat_data.values, aspect='auto', cmap=cmap,
                    norm=norm, interpolation='none')
ax_heat.set_yticks(range(len(classes_in_data)))
ax_heat.set_yticklabels([CLASS_LABELS[c] for c in classes_in_data], fontsize=6)
ax_heat.set_xticks([])
ax_heat.set_xlabel('Patients (sorted by P. aeruginosa abundance →)', fontsize=6, labelpad=4)

cbar = fig.colorbar(im, ax=ax_heat, fraction=0.028, pad=0.02, aspect=18)
cbar.set_label('Unique AROs', fontsize=5.5)
cbar.ax.tick_params(labelsize=5)

# ─── Lollipop — top 20 AROs ────────────────────────────────────────────────
top20 = prev.head(20).copy()
top20['color'] = top20['broad_class'].map(CLASS_COLORS).fillna('#BDBDBD')

def short_label(s):
    s = str(s)
    s = s.replace('Pseudomonas aeruginosa ', 'PA ')
    s = s.replace('Neisseria gonorrhoeae 23S rRNA with mutation conferring resistance to azithromycin',
                  'Ng 23S rRNA (azithromycin)')
    return s[:30]

top20['short'] = top20['gene_label'].apply(short_label)
top20 = top20.iloc[::-1].reset_index(drop=True)
y = np.arange(len(top20))

ax_lol.hlines(y, 0, top20['n_samples'], colors='#DDDDDD', linewidth=0.8, zorder=1)
ax_lol.scatter(top20['n_samples'], y, c=top20['color'], s=18, zorder=3, linewidths=0)
ax_lol.axvline(40, color='#AAAAAA', linewidth=0.5, linestyle='--')
ax_lol.set_yticks(y)
ax_lol.set_yticklabels(top20['short'], fontsize=5, fontstyle='italic')
ax_lol.set_xlabel('No. patients (of 40)', fontsize=7)
ax_lol.set_title('Most prevalent ARGs', fontsize=7, pad=4)
ax_lol.set_xlim(0, 44)
ax_lol.spines['top'].set_visible(False)
ax_lol.spines['right'].set_visible(False)
ax_lol.tick_params(axis='x', length=2.5, width=0.6, labelsize=6)
ax_lol.tick_params(axis='y', length=0)

leg_classes = sorted(top20['broad_class'].unique())
leg_h = [mpatches.Patch(color=CLASS_COLORS.get(c,'#BDBDBD'), label=CLASS_LABELS.get(c,c))
         for c in leg_classes]
ax_lol.legend(handles=leg_h, fontsize=5, frameon=False, loc='lower right',
              handlelength=0.9, labelspacing=0.25, ncol=1)

# ─── Mechanism donut ───────────────────────────────────────────────────────
# One row per unique ARO (each ARO has a single mechanism)
mech_counts = (df.drop_duplicates('ARO_id')['Resistance Mechanism']
               .fillna('unknown')
               .value_counts())
mech_simple = {}
for m, n in mech_counts.items():
    if m == 'antibiotic efflux':
        key = 'Efflux'
    elif 'inactivation' in m:
        key = 'Inactivation'
    elif 'alteration' in m and 'efflux' not in m:
        key = 'Target alteration'
    elif 'alteration' in m and 'efflux' in m:
        key = 'Target alt.+efflux'
    elif 'protection' in m:
        key = 'Target protection'
    elif 'replacement' in m:
        key = 'Target replacement'
    else:
        key = 'Other'
    mech_simple[key] = mech_simple.get(key, 0) + n

mech_df = pd.Series(mech_simple).sort_values(ascending=False)
colors_d = [MECH_COLORS.get(k, '#BDBDBD') for k in mech_df.index]
total_mech = mech_df.sum()

wedges, _, autotexts = ax_mech.pie(
    mech_df.values, labels=None, colors=colors_d,
    autopct=lambda p: f'{p:.0f}%' if p > 5 else '',
    startangle=90, pctdistance=0.78,
    wedgeprops={'linewidth': 0.5, 'edgecolor': 'white'},
)
for at in autotexts:
    at.set_fontsize(5.5); at.set_color('white'); at.set_fontweight('bold')

ax_mech.add_patch(plt.Circle((0, 0), 0.55, fc='white'))
ax_mech.text(0, 0, f'n={total_mech}\nAROs', ha='center', va='center',
             fontsize=6, fontweight='bold')
ax_mech.set_title('Resistance\nmechanism', fontsize=7, pad=4)

leg_mech = [mpatches.Patch(color=MECH_COLORS.get(k,'#BDBDBD'), label=k) for k in mech_df.index]
ax_mech.legend(handles=leg_mech, fontsize=5, frameon=False, loc='lower center',
               bbox_to_anchor=(0.5, -0.25), ncol=1, handlelength=0.9, labelspacing=0.2)

# ─── Stacked bar — group comparison ────────────────────────────────────────
pa_pos = [s for s in co_samples if meta.loc[s,'PA_sylph_pct'] > 0]
pa_neg = [s for s in co_samples if meta.loc[s,'PA_sylph_pct'] == 0]
azm_pos = [s for s in co_samples if meta.loc[s,'AZM'] == 'AZMpos']
azm_neg = [s for s in co_samples if meta.loc[s,'AZM'] == 'AZMneg']
leedsA  = [s for s in co_samples if meta.loc[s,'Leeds'] == 'LeedsA']
leedsB  = [s for s in co_samples if meta.loc[s,'Leeds'] == 'LeedsB']

groups = {
    f'All\n(n=40)':    co_samples,
    f'PA+\n(n={len(pa_pos)})':  pa_pos,
    f'PA−\n(n={len(pa_neg)})':  pa_neg,
    f'AZM+\n(n=20)':  azm_pos,
    f'AZM−\n(n=20)':  azm_neg,
}

x = np.arange(len(groups))
bottoms = np.zeros(len(groups))
for cls in classes_in_data:
    vals = np.array([
        pivot_co.loc[[s for s in samps if s in pivot_co.index], cls].mean()
        if any(s in pivot_co.index for s in samps) else 0
        for samps in groups.values()
    ])
    ax_bar.bar(x, vals, bottom=bottoms,
               color=CLASS_COLORS.get(cls, '#BDBDBD'), width=0.55, linewidth=0)
    bottoms += vals

# Separator lines
ax_bar.axvline(2.5, color='#AAAAAA', linewidth=0.8, linestyle=':')

# Significance annotations
bar_top = bottoms.max()
brac_y = bar_top * 1.03
ns_y   = bar_top * 1.06

# PA+ vs PA− (x=1, x=2): p < 0.0001 — significant
ax_bar.plot([1, 1, 2, 2], [brac_y, ns_y, ns_y, brac_y],
            color='#C0392B', linewidth=0.7, clip_on=False)
ax_bar.text(1.5, ns_y * 1.005, 'p < 0.0001', ha='center', va='bottom',
            fontsize=5.5, color='#C0392B')

# AZM+ vs AZM− (x=3, x=4): p = 0.239 — NS
ax_bar.plot([3, 3, 4, 4], [brac_y, ns_y, ns_y, brac_y],
            color='#555555', linewidth=0.7, clip_on=False)
ax_bar.text(3.5, ns_y * 1.005, 'p = 0.239', ha='center', va='bottom',
            fontsize=5.5, color='#888888')

ylim_top = bar_top * 1.18
ax_bar.set_ylim(0, ylim_top)

ax_bar.text(1.0, ylim_top * 0.98, 'PA pulmotype',
            ha='center', va='top', fontsize=6, color='#555')
ax_bar.text(3.5, ylim_top * 0.98, 'AZM status',
            ha='center', va='top', fontsize=6, color='#555')

ax_bar.set_xticks(x)
ax_bar.set_xticklabels(list(groups.keys()), fontsize=7)
ax_bar.set_ylabel('Mean unique AROs per patient', fontsize=7)
ax_bar.set_title('Resistome composition by clinical group', fontsize=7, pad=4)
ax_bar.spines['top'].set_visible(False)
ax_bar.spines['right'].set_visible(False)
ax_bar.tick_params(axis='both', length=2.5, width=0.6)

leg_all = [mpatches.Patch(color=CLASS_COLORS.get(c,'#BDBDBD'), label=CLASS_LABELS.get(c,c))
           for c in classes_in_data]
ax_bar.legend(handles=leg_all, fontsize=5.5, frameon=False, loc='upper right',
              handlelength=0.9, ncol=3, labelspacing=0.25, columnspacing=0.6)

# ─── Panel labels ──────────────────────────────────────────────────────────
ax_pa.text(-0.08, 1.05, 'A', transform=ax_pa.transAxes,
           fontsize=9, fontweight='bold', va='top')
ax_lol.text(-0.14, 1.02, 'B', transform=ax_lol.transAxes,
            fontsize=9, fontweight='bold', va='top')
ax_mech.text(-0.18, 1.02, 'C', transform=ax_mech.transAxes,
             fontsize=9, fontweight='bold', va='top')
ax_bar.text(-0.06, 1.04, 'D', transform=ax_bar.transAxes,
            fontsize=9, fontweight='bold', va='top')

# ─── Save ──────────────────────────────────────────────────────────────────
fig.savefig(OUTPDF, dpi=300)
fig.savefig(OUTPNG, dpi=300)
plt.close(fig)
print(f'Saved:\n  {OUTPDF}\n  {OUTPNG}')
