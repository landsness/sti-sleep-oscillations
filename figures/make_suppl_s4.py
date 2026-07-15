"""Supplemental Figure S5: Baseline Absolute SO vs Week 1 Behavior (null result)

Panel A: Baseline ROI1 Ipsi SO vs Week 1 Behavior
Panel B: Baseline ROI2 Ipsi SO vs Week 1 Behavior
Both are null results — baseline absolute SO does not predict recovery.
n=25, Spearman.
"""
import sys
from pathlib import Path
import matplotlib.patches as mpatches
import matplotlib.pyplot as plt
from scipy.stats import spearmanr

FIGURES_DIR = Path(__file__).parent
sys.path.insert(0, str(FIGURES_DIR))

from figure_style import (
    COLORS, FONT_SIZES, add_panel_label, format_stats,
    scatter_with_regression, save_figure, set_publication_style, load_data,
)

DATA_FILE = FIGURES_DIR.parent / 'WFCI AI Summary.xlsx'
OUTPUT_DIR = FIGURES_DIR / 'final'


def main():
    set_publication_style()
    data = load_data(DATA_FILE)
    bl, wk = data['bl'], data['wk']

    colors_all = [COLORS['sti_pos'] if bl.loc[m, 'secondary_thalamic_injury'] == 1
                  else COLORS['sti_neg'] for m in bl.index]
    patch_pos = mpatches.Patch(color=COLORS['sti_pos'], label='STI+')
    patch_neg = mpatches.Patch(color=COLORS['sti_neg'], label='STI-')

    y_wk = wk.loc[bl.index, 'behavior_score'].values

    panels = [
        ('so_power_roi1_ipsi', 'Baseline Ipsilateral SO Power\n' + r'(Lesion, $\times 10^{-5}$ µV²)', 'A'),
        ('so_power_roi2_ipsi', 'Baseline Ipsilateral SO Power\n' + r'(Perilesional, $\times 10^{-5}$ µV²)', 'B'),
    ]

    fig, axes = plt.subplots(1, 2, figsize=(10, 4.5))

    for ax, (col, xlabel, panel) in zip(axes, panels):
        x = bl[col].values * 1e5
        rho, p = spearmanr(x, y_wk)

        scatter_with_regression(
            ax, x, y_wk, colors_all,
            xlabel=xlabel,
            ylabel='Week 1 Forelimb Asymmetry Score (%)',
            stats_text=format_stats(rho=rho, p=p),
            loc='upper right',
        )
        ax.axhline(0, color='#666666', lw=0.8, ls='--', alpha=0.6, zorder=1)
        add_panel_label(ax, panel)

        if panel == 'A':
            ax.legend(handles=[patch_pos, patch_neg],
                      fontsize=FONT_SIZES['legend'], framealpha=0.8, loc='lower right')

    plt.tight_layout(pad=0.5)
    save_figure(fig, 'SupplFig_S4_BaselineSO_Null', OUTPUT_DIR)
    plt.close()
    print('SupplFig S4 done.')


if __name__ == '__main__':
    main()
