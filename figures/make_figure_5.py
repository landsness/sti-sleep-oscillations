"""make_figure_6.py -- Figure 5: Pre-Stroke Network State Predicts Recovery

Panel A: Baseline Laterality Index (ROI1) vs Week 1 Behavior  [n=25, Spearman]
Panel B: Baseline Laterality Index (ROI1) vs Infarct Volume   [n=25, Spearman]
         Independence check -- LI is on X in both panels for parallel structure
"""

import sys
from pathlib import Path

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from scipy.stats import spearmanr, rankdata, pearsonr
from scipy.stats import linregress

FIGURES_DIR = Path(__file__).parent
sys.path.insert(0, str(FIGURES_DIR))

from figure_style import (
    COLORS, FONT_SIZES, LINE_WIDTHS, MARKER_SIZES,
    set_publication_style, add_panel_label, add_stats_box,
    scatter_with_regression, format_pval, format_stats, save_figure,
    load_data,
)

DATA_FILE = FIGURES_DIR.parent / 'WFCI AI Summary.xlsx'
OUT_DIR   = FIGURES_DIR / 'final'


def load_figure_data():
    return load_data(DATA_FILE)


def partial_spear(x, y, z):
    """Partial Spearman correlation of x~y controlling for z (rank-residual)."""
    x, y, z = np.asarray(x, float), np.asarray(y, float), np.asarray(z, float)
    mask = ~(np.isnan(x) | np.isnan(y) | np.isnan(z))
    rx = rankdata(x[mask]).astype(float)
    ry = rankdata(y[mask]).astype(float)
    rz = rankdata(z[mask]).astype(float)
    def resid(a, b):
        sl, ic, *_ = linregress(b, a)
        return a - (sl * b + ic)
    return pearsonr(resid(rx, rz), resid(ry, rz))


def prepare_data(data):
    """Use full n=25 dataset -- outlier exclusion only applies to SO power analyses."""
    bl = data['bl']
    wk = data['wk']   # all 25 animals

    # Align on shared index (should be identical for bl and wk)
    shared = bl.index.intersection(wk.index)
    bl = bl.loc[shared]
    wk = wk.loc[shared]

    x_li     = bl['LI_roi1'].values.astype(float)
    x_stroke = bl['stroke_size'].values.astype(float)
    y_wk1    = wk['behavior_score'].values.astype(float)
    n        = len(shared)

    colors_all = [
        COLORS['sti_pos'] if bl['secondary_thalamic_injury'].iloc[i] == 1
        else COLORS['sti_neg']
        for i in range(n)
    ]

    return dict(x_li=x_li, x_stroke=x_stroke, y_wk1=y_wk1,
                n=n, point_colors=colors_all)


def main():
    set_publication_style()
    data = load_figure_data()
    pd_  = prepare_data(data)

    x_li, x_stroke, y_wk1 = pd_['x_li'], pd_['x_stroke'], pd_['y_wk1']
    n, point_colors = pd_['n'], pd_['point_colors']

    # Panel A: Spearman LI vs behavior + partial rho
    rho_a, p_a = spearmanr(x_li, y_wk1)
    partial_rho_a, partial_p_a = partial_spear(x_li, y_wk1, x_stroke)

    stats_text_a = (
        f'ρ = {rho_a:.3f}\n'          # ρ
        f'{format_pval(p_a)}\n'
        f'partial ρ = {partial_rho_a:.3f}\n'
        f'{format_pval(partial_p_a)}'
    )

    # Panel B: Spearman LI vs infarct (independence check)
    # X = LI (same as Panel A), Y = infarct volume
    rho_b, p_b = spearmanr(x_li, x_stroke)

    stats_text_b = (
        f'ρ = {rho_b:.3f}\n'
        f'{format_pval(p_b)}'
    )

    patch_pos = mpatches.Patch(color=COLORS['sti_pos'], label='STI+')
    patch_neg = mpatches.Patch(color=COLORS['sti_neg'], label='STI−')

    fig, axes = plt.subplots(1, 2, figsize=(10, 4.8))

    # Panel A: LI (x) vs Week 1 behavior (y)
    scatter_with_regression(
        axes[0],
        x=x_li, y=y_wk1,
        colors=point_colors,
        xlabel='Baseline Laterality Index (Lesion)',
        ylabel='Week 1 Forelimb Asymmetry Score (%)',
        stats_text=stats_text_a,
        loc='lower left',
    )
    axes[0].axhline(0, color='#666666', lw=0.8, ls='--', alpha=0.6, zorder=1)
    axes[0].legend(handles=[patch_pos, patch_neg],
                   fontsize=FONT_SIZES['legend'], framealpha=0.85, loc='upper right')
    add_panel_label(axes[0], 'A')

    # Panel B: LI (x) vs Infarct Volume (y) — both axes now parallel to Panel A
    scatter_with_regression(
        axes[1],
        x=x_li, y=x_stroke,
        colors=point_colors,
        xlabel='Baseline Laterality Index (Lesion)',
        ylabel='Infarct Volume (mm³)',   # mm³
        stats_text=stats_text_b,
        loc='upper right',
    )
    add_panel_label(axes[1], 'B')

    plt.tight_layout(pad=0.5)
    save_figure(fig, 'Figure5_BaselineLI', OUT_DIR)
    plt.close()
    print(f'Figure 5 done. n={n}, rho_A={rho_a:.3f} p={p_a:.3f}, '
          f'partial_rho={partial_rho_a:.3f} p={partial_p_a:.3f}')


if __name__ == '__main__':
    main()
