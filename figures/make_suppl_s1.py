"""Supplemental Figure S1: Infarct Size vs Behavior (standalone clean script)."""
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
    bl, ac, wk = data['bl'], data['ac'], data['wk']

    colors_all = [COLORS['sti_pos'] if bl.loc[m, 'secondary_thalamic_injury'] == 1
                  else COLORS['sti_neg'] for m in bl.index]
    patch_pos = mpatches.Patch(color=COLORS['sti_pos'], label='STI+')
    patch_neg = mpatches.Patch(color=COLORS['sti_neg'], label='STI-')

    x = bl['stroke_size'].values
    y_ac = ac.loc[bl.index, 'behavior_score'].values
    y_wk = wk.loc[bl.index, 'behavior_score'].values

    r_ac, p_ac = spearmanr(x, y_ac)
    r_wk, p_wk = spearmanr(x, y_wk)

    fig, axes = plt.subplots(1, 2, figsize=(10, 4.5))

    # Panel A: Infarct vs acute behavior
    scatter_with_regression(
        axes[0], x, y_ac, colors_all,
        xlabel='Infarct Volume (mm³)',
        ylabel='Acute Forelimb Asymmetry Score (%)',
        stats_text=format_stats(rho=r_ac, p=p_ac),
        loc='upper left',
    )
    axes[0].axhline(0, color='#666666', lw=0.8, ls='--', alpha=0.6, zorder=1)
    axes[0].legend(handles=[patch_pos, patch_neg],
                   fontsize=FONT_SIZES['legend'], framealpha=0.8, loc='lower right')
    add_panel_label(axes[0], 'A')

    # Panel B: Infarct vs week 1 behavior — stats lower left (data trends upper right)
    scatter_with_regression(
        axes[1], x, y_wk, colors_all,
        xlabel='Infarct Volume (mm³)',
        ylabel='Week 1 Forelimb Asymmetry Score (%)',
        stats_text=format_stats(rho=r_wk, p=p_wk),
        loc='lower right',
    )
    axes[1].axhline(0, color='#666666', lw=0.8, ls='--', alpha=0.6, zorder=1)
    add_panel_label(axes[1], 'B')

    plt.tight_layout(pad=0.5, w_pad=2.0)
    save_figure(fig, 'SupplFig_S1_InfarctVsBehavior', OUTPUT_DIR)
    plt.close()
    print('SupplFig S1 done.')


if __name__ == '__main__':
    main()
