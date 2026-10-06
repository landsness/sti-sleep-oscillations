"""Supplementary Figure 7: Baseline Absolute SO vs Week 1 Behavior (null result)

2 x 2, one panel per region x hemisphere (rev. 2026-10-05: contralateral panels C, D
added so every comparison cited to this figure in Results 3.4 is plotted):
  A lesion ipsilateral     B perilesional ipsilateral
  C lesion contralateral   D perilesional contralateral
All are null results: absolute baseline SO power does not predict week-1 outcome.
n = 25, Spearman rho with uncorrected p (exploratory family; q in Suppl Table 2).
"""
import sys
from pathlib import Path
import matplotlib.patches as mpatches
import matplotlib.pyplot as plt
from scipy.stats import spearmanr

FIGURES_DIR = Path(__file__).resolve().parent
sys.path[:0] = [str(FIGURES_DIR), str(FIGURES_DIR.parent), str(FIGURES_DIR.parent / 'analysis')]
import config  # noqa: E402

from figure_style import (
    COLORS, FONT_SIZES, add_panel_label, format_stats,
    scatter_with_regression, save_figure, set_publication_style, load_data,
)

DATA_FILE = config.XLSX
OUTPUT_DIR = config.FIG_OUT


def main():
    set_publication_style()
    data = load_data(DATA_FILE)
    bl, wk = data['bl'], data['wk']

    colors_all = [COLORS['sti_pos'] if bl.loc[m, 'secondary_thalamic_injury'] == 1
                  else COLORS['sti_neg'] for m in bl.index]
    patch_pos = mpatches.Patch(color=COLORS['sti_pos'], label='STI+')
    patch_neg = mpatches.Patch(color=COLORS['sti_neg'], label='STI−')

    y_wk = wk.loc[bl.index, 'behavior_score'].values

    unit = r'$\times 10^{-5}$ (ΔF/F)²/Hz'
    panels = [
        ('so_power_roi1_ipsi',   f'Baseline Ipsilateral SO Power\n(Lesion, {unit})',         'A'),
        ('so_power_roi2_ipsi',   f'Baseline Ipsilateral SO Power\n(Perilesional, {unit})',   'B'),
        ('so_power_roi1_contra', f'Baseline Contralateral SO Power\n(Lesion, {unit})',       'C'),
        ('so_power_roi2_contra', f'Baseline Contralateral SO Power\n(Perilesional, {unit})', 'D'),
    ]

    fig, axes = plt.subplots(2, 2, figsize=(9.2, 8.4), sharey=True)

    for ax, (col, xlabel, panel) in zip(axes.ravel(), panels):
        x = bl[col].values * 1e5
        rho, p = spearmanr(x, y_wk)
        print(f'{col:22s} rho={rho:+.3f} p={p:.3f}')

        scatter_with_regression(
            ax, x, y_wk, colors_all,
            xlabel=xlabel,
            ylabel='Week 1 Forelimb Asymmetry Score (%)' if panel in 'AC' else '',
            stats_text=format_stats(rho=rho, p=p),
            loc='upper left' if panel == 'A' else 'upper right',   # A: keep the text off the (20.8, 55) point
        )
        ax.axhline(0, color='#666666', lw=0.8, ls='--', alpha=0.6, zorder=1)
        ax.set_xlim(left=0)
        add_panel_label(ax, panel)

        if panel == 'A':
            ax.legend(handles=[patch_pos, patch_neg],
                      fontsize=FONT_SIZES['legend'], framealpha=0.8, loc='lower right')

    plt.tight_layout(pad=0.5)
    save_figure(fig, 'SupplFig7_BaselineSO_Null', OUTPUT_DIR)
    plt.close()
    print('Supplementary Figure 7 done.')


if __name__ == '__main__':
    main()
