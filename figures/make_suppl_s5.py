"""Supplemental Figure S5: Baseline LI ROI2 Null (single-panel scatter)

Panel A: Baseline Laterality Index ROI2 vs Week 1 Behavior (null result)
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
    COLORS, FONT_SIZES, add_panel_label, format_pval,
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

    x_li_r2 = bl['LI_roi2'].values
    y_wk = wk.loc[bl.index, 'behavior_score'].values

    rho, p = spearmanr(x_li_r2, y_wk)
    stats_text = f'ρ = {rho:.3f}\n{format_pval(p)}'

    fig, ax = plt.subplots(1, 1, figsize=(5, 4.8))

    scatter_with_regression(
        ax, x_li_r2, y_wk, colors_all,
        xlabel='Baseline Laterality Index (Perilesional)',
        ylabel='Week 1 Forelimb Asymmetry Score (%)',
        stats_text=stats_text,
        loc='lower right',
    )
    ax.axhline(0, color='#666666', lw=0.8, ls='--', alpha=0.6, zorder=1)
    ax.legend(handles=[patch_pos, patch_neg],
              fontsize=FONT_SIZES['legend'], framealpha=0.8, loc='upper left')
    add_panel_label(ax, 'A')

    plt.tight_layout(pad=0.5)
    save_figure(fig, 'SupplFig_S5_BaselineLI_ROI2_Null', OUTPUT_DIR)
    plt.close()
    print('SupplFig S5 done.')


if __name__ == '__main__':
    main()
