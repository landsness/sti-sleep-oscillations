"""R1.13 panels: infarct volume x acute (24h) SO power, STI-colored, mock style.
Four scatters: lesion ipsi (ROI1), peri ipsi (ROI2), lesion contra (ROI1), peri contra (ROI2).
Uses the manuscript's figure_style for exact visual consistency with Figure 2.
"""
import sys
from pathlib import Path
import numpy as np, pandas as pd
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from scipy.stats import spearmanr

FS = '/mnt/user-data/uploads/STI_SO/STI_SO_figures/figures'
sys.path.insert(0, FS)
from figure_style import (COLORS, FONT_SIZES, set_publication_style,
                          scatter_with_regression, add_panel_label, format_pval)

DATA = '/mnt/user-data/uploads/STI_SO/data/WFCI AI Summary.xlsx'
OUT = '/mnt/user-data/outputs'

df = pd.read_excel(DATA, sheet_name='AI Summary '); df.columns = df.columns.str.strip()
for c in ['mouse_id', 'time_point']: df[c] = df[c].astype(str).str.strip()
bl = df[df.time_point == 'baseline'].set_index('mouse_id').copy()
ac = df[df.time_point == '24 hours'].set_index('mouse_id').copy().loc[bl.index]

infarct = bl['stroke_size'].values.astype(float)
sti = bl['secondary_thalamic_injury'].values
pt_colors = [COLORS['sti_pos'] if s == 1 else COLORS['sti_neg'] for s in sti]
SCALE = 1e-5  # display SO power in x10^-5 units, matching the mock

PANELS = [
    ('so_power_roi1_ipsi',   'Acute Ipsilateral SO Power\n(ROI1, ×10$^{-5}$ µV²)',   'lesion ipsi'),
    ('so_power_roi2_ipsi',   'Acute Ipsilateral SO Power\n(ROI2, ×10$^{-5}$ µV²)',   'peri ipsi'),
    ('so_power_roi1_contra', 'Acute Contralateral SO Power\n(ROI1, ×10$^{-5}$ µV²)', 'lesion contra'),
    ('so_power_roi2_contra', 'Acute Contralateral SO Power\n(ROI2, ×10$^{-5}$ µV²)', 'peri contra'),
]

def stats_text(col):
    y = ac[col].values.astype(float)
    rho, p = spearmanr(infarct, y)
    return f'ρ = {rho:.3f}\n{format_pval(p)}', rho, p

set_publication_style()
patch_pos = mpatches.Patch(color=COLORS['sti_pos'], label='STI+')
patch_neg = mpatches.Patch(color=COLORS['sti_neg'], label='STI−')

# ---- Combined 2x2 review figure ----
fig, axes = plt.subplots(2, 2, figsize=(9.2, 8.4))
labels = ['A', 'B', 'C', 'D']
for ax, (col, ylab, _), lab in zip(axes.ravel(), PANELS, labels):
    txt, rho, p = stats_text(col)
    y = ac[col].values.astype(float) / SCALE
    loc = 'upper right' if rho <= 0 else 'lower right'
    scatter_with_regression(ax, x=infarct, y=y, colors=pt_colors,
                            xlabel='Infarct Volume (mm³)', ylabel=ylab,
                            stats_text=txt, loc=loc)
    add_panel_label(ax, lab)
axes.ravel()[0].legend(handles=[patch_pos, patch_neg], fontsize=FONT_SIZES['legend'],
                       framealpha=0.85, loc='lower left')
plt.tight_layout(pad=0.6)
fig.savefig(f'{OUT}/R1.13_infarct_vs_acuteSO_2x2.png', dpi=200, bbox_inches='tight')
plt.close()

# ---- Individual panels (for slotting into Figure 2 / supplement) ----
for col, ylab, tag in PANELS:
    txt, rho, p = stats_text(col)
    y = ac[col].values.astype(float) / SCALE
    loc = 'upper right' if rho <= 0 else 'lower right'
    fig, ax = plt.subplots(figsize=(4.4, 4.2))
    scatter_with_regression(ax, x=infarct, y=y, colors=pt_colors,
                            xlabel='Infarct Volume (mm³)', ylabel=ylab,
                            stats_text=txt, loc=loc)
    if tag == 'lesion ipsi':
        ax.legend(handles=[patch_pos, patch_neg], fontsize=FONT_SIZES['legend'],
                  framealpha=0.85, loc='lower left')
    plt.tight_layout(pad=0.5)
    slug = tag.replace(' ', '_')
    fig.savefig(f'{OUT}/R1.13_infarct_{slug}.png', dpi=300, bbox_inches='tight')
    plt.close()
    print(f'{tag:16s} rho={rho:+.3f} p={p:.4f}')
print('done')
