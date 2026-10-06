"""Supplementary Figure 3 — Infarct volume vs. acute (24 h) slow-oscillation power.

Reviewer 1.13. Four STI-colored scatterplots (infarct volume on x):
  A  lesional ipsilateral   (ROI1)   rho = -0.809, p < 0.001
  B  perilesional ipsilateral (ROI2) rho = -0.742, p < 0.001
  C  lesional contralateral  (ROI1)  rho = +0.092, p = 0.663  (null / specificity)
  D  perilesional contralateral (ROI2) rho = -0.078, p = 0.712 (null / specificity)

Acute ipsilateral SO power scales with lesion size, ipsilaterally and
region-wide, while contralateral SO does not — the continuous-severity
counterpart to the STI (categorical) result in Figure 2E, and support for
lesional-ROI reliability (larger lesion -> more suppression, not a saturating
artifact). Uses figure_style for exact visual consistency with the main figures.

Output: results/figures/{svg,pdf,png}/SupplFig3_InfarctVsAcuteSO.*
"""
import sys
from pathlib import Path
import numpy as np, pandas as pd
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from scipy.stats import spearmanr

FIGURES_DIR = Path(__file__).resolve().parent
sys.path[:0] = [str(FIGURES_DIR), str(FIGURES_DIR.parent), str(FIGURES_DIR.parent / 'analysis')]
import config  # noqa: E402
from figure_style import (COLORS, FONT_SIZES, set_publication_style,
                          scatter_with_regression, add_panel_label, format_pval,
                          save_figure)

from figure_style import load_data, q_value, format_qval  # noqa: E402

SCALE = 1e-5  # display SO power in x10^-5 units, matching Figure 2

PANELS = [
    ('so_power_roi1_ipsi',   'Acute Ipsilateral SO Power\n(Lesion, ×10$^{-5}$ (ΔF/F)²/Hz)',   'lesion ipsi'),
    ('so_power_roi2_ipsi',   'Acute Ipsilateral SO Power\n(Perilesional, ×10$^{-5}$ (ΔF/F)²/Hz)',   'peri ipsi'),
    ('so_power_roi1_contra', 'Acute Contralateral SO Power\n(Lesion, ×10$^{-5}$ (ΔF/F)²/Hz)', 'lesion contra'),
    ('so_power_roi2_contra', 'Acute Contralateral SO Power\n(Perilesional, ×10$^{-5}$ (ΔF/F)²/Hz)', 'peri contra'),
]


def main():
    data = load_data(config.XLSX)
    bl = data['bl']
    ac = data['ac'].loc[bl.index]
    infarct = bl['stroke_size'].values.astype(float)
    sti = bl['secondary_thalamic_injury'].values
    pt_colors = [COLORS['sti_pos'] if s == 1 else COLORS['sti_neg'] for s in sti]

    set_publication_style()
    patch_pos = mpatches.Patch(color=COLORS['sti_pos'], label='STI+')
    patch_neg = mpatches.Patch(color=COLORS['sti_neg'], label='STI−')

    fig, axes = plt.subplots(2, 2, figsize=(9.2, 8.4))
    for ax, (col, ylab, tag), lab in zip(axes.ravel(), PANELS, ['A', 'B', 'C', 'D']):
        rho, p = spearmanr(infarct, ac[col].values.astype(float))
        q = q_value(f'infarct_vs_acute_so_{tag}', p)
        y = ac[col].values.astype(float) / SCALE
        scatter_with_regression(ax, x=infarct, y=y, colors=pt_colors,
                                xlabel='Infarct Volume (mm³)', ylabel=ylab,
                                stats_text=f'ρ = {rho:.3f}\n{format_pval(p)}\n{format_qval(q)}'.replace('-', '\u2212'),
                                loc='upper right' if rho <= 0 else 'lower right')
        ax.set_ylim(bottom=0)  # SO power is non-negative
        add_panel_label(ax, lab)
        print(f'{tag:16s} rho={rho:+.3f} p={p:.4f}')
    axes.ravel()[0].legend(handles=[patch_pos, patch_neg], fontsize=FONT_SIZES['legend'],
                           framealpha=0.85, loc='center right')
    plt.tight_layout(pad=0.6)
    save_figure(fig, 'SupplFig3_InfarctVsAcuteSO', config.FIG_OUT)
    plt.close(fig)


if __name__ == '__main__':
    main()
