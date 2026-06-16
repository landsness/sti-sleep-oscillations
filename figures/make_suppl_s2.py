"""Supplemental Figure S2: STI-Subgroup Robustness of the SO-Behavior Relationship

Panel A: STI+ subgroup — Acute ROI1 Ipsi SO vs Acute Behavior
Panel B: STI- subgroup — Acute ROI1 Ipsi SO vs Acute Behavior

Shows that the SO-behavior relationship holds within each STI group separately.
Content from old make_suppl_s4.py.
"""
import sys
from pathlib import Path
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
    bl, ac = data['bl'], data['ac']
    sti_pos, sti_neg = data['sti_pos'], data['sti_neg']

    fig, axes = plt.subplots(1, 2, figsize=(10, 4.5))

    for ax, mask, color, panel in [
        (axes[0], sti_pos, COLORS['sti_pos'], 'A'),
        (axes[1], sti_neg, COLORS['sti_neg'], 'B'),
    ]:
        sub_ac = ac[mask]
        x = sub_ac['so_power_roi1_ipsi'].values * 1e5
        y = sub_ac['behavior_score'].values
        rho, p = spearmanr(x, y)

        scatter_with_regression(
            ax, x, y, [color] * len(x),
            xlabel=r'Acute Ipsilateral SO Power (ROI1, $\times 10^{-5}$ µV²)',
            ylabel='Acute Forelimb Asymmetry Score (%)',
            stats_text=format_stats(rho=rho, p=p),
            loc='upper right',
        )
        ax.axhline(0, color='#666666', lw=0.8, ls='--', alpha=0.6, zorder=1)
        add_panel_label(ax, panel)

    plt.tight_layout(pad=0.5)
    save_figure(fig, 'SupplFig_S2_STISubgroups', OUTPUT_DIR)
    plt.close()
    print('SupplFig S2 done.')


if __name__ == '__main__':
    main()
