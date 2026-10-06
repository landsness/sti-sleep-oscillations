"""supplfig6_wk1_so_vs_behavior.py -- Supplementary Figure 6 (Reviewer 1.16; replaces Table 1).

Week-1 SO power (% of each animal's baseline) vs week-1 forelimb asymmetry,
one panel per region x hemisphere, n = 23 (the two week-1-excluded animals
omitted; see Methods). Points coloured by STI group; a single pooled
regression line per panel. Stats: Spearman rho, uncorrected p and BH q from
the confirmatory correlation family (analysis/fdr.py).

  A lesion ipsilateral     B perilesional ipsilateral
  C lesion contralateral   D perilesional contralateral

Output: results/figures/{png,pdf,svg}/SupplFig6_Wk1SO_vs_Behavior.*
"""
import sys
from pathlib import Path

import matplotlib.patches as mpatches
import matplotlib.pyplot as plt
from scipy.stats import spearmanr

FIGURES_DIR = Path(__file__).resolve().parent
sys.path[:0] = [str(FIGURES_DIR), str(FIGURES_DIR.parent), str(FIGURES_DIR.parent / 'analysis')]
import config  # noqa: E402

from figure_style import (  # noqa: E402
    COLORS, FONT_SIZES, add_panel_label, format_stats, q_value,
    scatter_with_regression, save_figure, set_publication_style, load_data,
)

PANELS = [
    ('lesion ipsi',   'so_power_roi1_ipsi',   'Week 1 Ipsilateral SO Power\n(Lesion, % Baseline)',          'A'),
    ('peri ipsi',     'so_power_roi2_ipsi',   'Week 1 Ipsilateral SO Power\n(Perilesional, % Baseline)',    'B'),
    ('lesion contra', 'so_power_roi1_contra', 'Week 1 Contralateral SO Power\n(Lesion, % Baseline)',        'C'),
    ('peri contra',   'so_power_roi2_contra', 'Week 1 Contralateral SO Power\n(Perilesional, % Baseline)',  'D'),
]


def main():
    set_publication_style()
    data = load_data(config.XLSX)
    bl, wkc = data['bl'], data['wk_clean']
    idx = wkc.index                                   # n = 23
    y = wkc['behavior_score'].values
    colors = [COLORS['sti_pos'] if bl.loc[m, 'secondary_thalamic_injury'] == 1 else COLORS['sti_neg']
              for m in idx]
    patch_pos = mpatches.Patch(color=COLORS['sti_pos'], label='STI+')
    patch_neg = mpatches.Patch(color=COLORS['sti_neg'], label='STI−')

    fig, axes = plt.subplots(2, 2, figsize=(9.2, 8.4))
    for ax, (key, col, xlabel, lab) in zip(axes.ravel(), PANELS):
        x = (wkc[col] / bl.loc[idx, col] * 100).values
        rho, p = spearmanr(x, y)
        q = q_value(f'wk1_so_vs_wk1_behavior_{key}', p)
        scatter_with_regression(ax, x, y, colors, xlabel=xlabel,
                                ylabel='Week 1 Forelimb Asymmetry Score (%)',
                                stats_text=format_stats(rho=rho, p=p, q=q), loc='upper right')
        ax.axhline(0, color='#666666', lw=0.8, ls='--', alpha=0.6, zorder=1)
        ax.axvline(100, color='#999999', lw=0.8, ls=':', alpha=0.7, zorder=1)   # baseline = 100%
        add_panel_label(ax, lab)
        print(f'{key:14s} rho={rho:+.3f} p={p:.3f} q={q:.3f}')
    axes[0, 0].legend(handles=[patch_pos, patch_neg], fontsize=FONT_SIZES['legend'],
                      framealpha=0.85, loc='lower right')
    plt.tight_layout(pad=0.6)
    save_figure(fig, 'SupplFig6_Wk1SO_vs_Behavior', config.FIG_OUT)
    plt.close(fig)


if __name__ == '__main__':
    main()
