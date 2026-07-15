"""make_figure_2.py -- Figure 2: Acute SO Suppression (7-panel, 3-row)

Row 0 — Overview (strip+box):
  A: Lateralization + spatial gradient (all animals, 4 channels)
  B: STI group difference in ipsilateral SO (4 channels)

Row 1 — Ipsilateral SO predicts behavior (scatter):
  C: ROI1 ipsi vs acute behavior
  D: ROI2 ipsi vs acute behavior
  E: Acute LI (ROI1) vs acute behavior

Row 2 — Contralateral SO null (scatter):
  F: ROI1 contra vs acute behavior
  G: ROI2 contra vs acute behavior
"""

import sys
from pathlib import Path

import matplotlib.patches as mpatches
import matplotlib.pyplot as plt
import numpy as np
from scipy import stats as scipy_stats
from scipy.stats import mannwhitneyu, wilcoxon

FIGURES_DIR = Path(__file__).parent
sys.path.insert(0, str(FIGURES_DIR))

from figure_style import (
    COLORS, FONT_SIZES, LINE_WIDTHS, MARKER_SIZES,
    add_panel_label, format_stats, format_pval, load_data,
    save_figure, scatter_with_regression, set_publication_style,
    style_ax,
)

DATA_FILE = FIGURES_DIR.parent / 'WFCI AI Summary.xlsx'
OUTPUT_DIR = FIGURES_DIR / 'final'

# Axis labels
SO_IPSI_ROI1  = r'Acute Ipsilateral SO Power' + '\n' + r'(Lesion, $\times 10^{-5}$ µV²)'
SO_IPSI_ROI2  = r'Acute Ipsilateral SO Power' + '\n' + r'(Perilesional, $\times 10^{-5}$ µV²)'
SO_CONTRA_ROI1 = r'Acute Contralateral SO Power' + '\n' + r'(Lesion, $\times 10^{-5}$ µV²)'
SO_CONTRA_ROI2 = r'Acute Contralateral SO Power' + '\n' + r'(Perilesional, $\times 10^{-5}$ µV²)'
LI_ROI1_LABEL = 'Acute Laterality Index (ROI1)'
BEH_ACUTE = 'Acute Forelimb\nAsymmetry Score (%)'

# Channel order: ipsilateral left, contralateral right
COLS = ['so_power_roi1_ipsi', 'so_power_roi2_ipsi',
        'so_power_roi1_contra', 'so_power_roi2_contra']
CH_LABELS = ['Lesion\nIpsi', 'Peri\nIpsi', 'Lesion\nContra', 'Peri\nContra']
CH_COLORS  = [COLORS['roi1'], COLORS['roi2'], COLORS['roi1'], COLORS['roi2']]


def spearman_stats(x, y):
    mask = np.isfinite(np.asarray(x, float)) & np.isfinite(np.asarray(y, float))
    xm, ym = np.asarray(x, float)[mask], np.asarray(y, float)[mask]
    rho, p = scipy_stats.spearmanr(xm, ym)
    return rho, p


def strip_box(ax, xi, vals, color, width=0.28, jitter_scale=0.06):
    np.random.seed(int(xi * 100) % 9999)
    bp = ax.boxplot([vals], positions=[xi], widths=width,
                    patch_artist=True,
                    medianprops=dict(color='black', lw=LINE_WIDTHS['median']),
                    whiskerprops=dict(lw=LINE_WIDTHS['box_whisker']),
                    capprops=dict(lw=LINE_WIDTHS['box_whisker']),
                    boxprops=dict(lw=LINE_WIDTHS['box_whisker']),
                    flierprops=dict(marker=''), zorder=2)
    bp['boxes'][0].set_facecolor(color)
    bp['boxes'][0].set_alpha(0.25)
    jit = np.random.uniform(-jitter_scale, jitter_scale, len(vals))
    ax.scatter(xi + jit, vals, color=color, s=MARKER_SIZES['strip'],
               alpha=0.75, edgecolors='white', lw=0.3, zorder=3)


def sig_bracket(ax, x1, x2, y, p, dy=3):
    sig = p < 0.05
    color = 'black' if sig else '#999999'
    fw = 'bold' if sig else 'normal'
    ax.plot([x1, x2], [y, y], color=color, lw=1.0, zorder=5)
    ax.plot([x1, x1], [y - dy*0.4, y], color=color, lw=1.0, zorder=5)
    ax.plot([x2, x2], [y - dy*0.4, y], color=color, lw=1.0, zorder=5)
    ax.text((x1+x2)/2, y + dy*0.3, format_pval(p),
            ha='center', va='bottom', fontsize=FONT_SIZES['stats'],
            fontweight=fw, color=color)


def main():
    set_publication_style()
    data = load_data(DATA_FILE)
    bl, ac = data['bl'], data['ac']
    sti_pos, sti_neg = data['sti_pos'], data['sti_neg']
    sti_pos_idx = bl.index[sti_pos]
    sti_neg_idx = bl.index[sti_neg]

    colors_all = [COLORS['sti_pos'] if sti_pos.get(mid, False) else COLORS['sti_neg']
                  for mid in ac.index]

    y_ac = ac['behavior_score'].values
    pct = {c: (ac[c] / bl[c] * 100).values for c in COLS}

    patch_pos = mpatches.Patch(color=COLORS['sti_pos'], label='STI+')
    patch_neg = mpatches.Patch(color=COLORS['sti_neg'], label='STI-')
    patch_r1  = mpatches.Patch(color=COLORS['roi1'], label='Lesion')
    patch_r2  = mpatches.Patch(color=COLORS['roi2'], label='Perilesional')

    # 2×3 grid — 6 panels:
    # Row 0: A (strip+box all animals), B (strip+box STI groups), C (ROI1 ipsi scatter)
    # Row 1: D (ROI2 ipsi scatter),     E (ROI1 contra scatter),   F (ROI2 contra scatter)
    fig, axes = plt.subplots(2, 3, figsize=(15, 9))

    # ── Panel A: All animals, 4 channels — lateralization + spatial gradient ──────
    ax = axes[0, 0]
    x_pos = np.array([0, 1, 2.4, 3.4])
    for xi, col, color in zip(x_pos, COLS, CH_COLORS):
        strip_box(ax, xi, pct[col], color)

    ax.axhline(100, color='#888888', lw=0.9, ls='--', alpha=0.6)
    v = [pct[c] for c in COLS]
    _, p_r1_hemi  = wilcoxon(pct['so_power_roi1_ipsi'], pct['so_power_roi1_contra'])
    _, p_r2_hemi  = wilcoxon(pct['so_power_roi2_ipsi'], pct['so_power_roi2_contra'])
    _, p_gradient = wilcoxon(pct['so_power_roi1_ipsi'], pct['so_power_roi2_ipsi'])
    ymax = max(vv.max() for vv in v)
    sig_bracket(ax, x_pos[0], x_pos[1], ymax + 8,  p_gradient, dy=3)
    sig_bracket(ax, x_pos[0], x_pos[2], ymax + 28, p_r1_hemi,  dy=3)
    sig_bracket(ax, x_pos[1], x_pos[3], ymax + 50, p_r2_hemi,  dy=3)
    ax.set_xticks(x_pos)
    ax.set_xticklabels(CH_LABELS, fontsize=FONT_SIZES['tick'])
    ax.set_ylabel('SO Power (% Baseline)', fontsize=FONT_SIZES['axis_label'])
    ax.set_ylim(-10, ymax + 90)
    style_ax(ax)
    add_panel_label(ax, 'A')

    # ── Panel B: STI+ vs STI- for all 4 channels ─────────────────────────────────
    ax = axes[0, 1]
    pair_centers = np.array([0, 1, 2.4, 3.4])
    w = 0.28
    for xi, col in zip(pair_centers, COLS):
        vpos = (ac.loc[sti_pos_idx, col] / bl.loc[sti_pos_idx, col] * 100).values
        vneg = (ac.loc[sti_neg_idx, col] / bl.loc[sti_neg_idx, col] * 100).values
        strip_box(ax, xi - w*0.65, vpos, COLORS['sti_pos'], width=w)
        strip_box(ax, xi + w*0.65, vneg, COLORS['sti_neg'], width=w)
        _, p_mw = mannwhitneyu(vpos, vneg, alternative='two-sided')
        y_top = max(vpos.max(), vneg.max()) + 5
        sig_bracket(ax, xi - w*0.65, xi + w*0.65, y_top, p_mw, dy=3)

    ax.axhline(100, color='#888888', lw=0.9, ls='--', alpha=0.6)
    ax.set_xticks(pair_centers)
    ax.set_xticklabels(CH_LABELS, fontsize=FONT_SIZES['tick'])
    ax.set_ylabel('SO Power (% Baseline)', fontsize=FONT_SIZES['axis_label'])
    ax.legend(handles=[patch_pos, patch_neg], fontsize=FONT_SIZES['legend'],
              loc='upper left', framealpha=0.9)
    style_ax(ax)
    add_panel_label(ax, 'B')

    # ── Panel C: ROI1 ipsi vs acute behavior ──────────────────────────────────────
    rho, p = spearman_stats(ac['so_power_roi1_ipsi'].values * 1e5, y_ac)
    scatter_with_regression(
        axes[0, 2], ac['so_power_roi1_ipsi'].values * 1e5, y_ac, colors_all,
        xlabel=SO_IPSI_ROI1, ylabel=BEH_ACUTE,
        stats_text=format_stats(rho=rho, p=p), loc='upper right',
    )
    axes[0, 2].legend(handles=[patch_pos, patch_neg], fontsize=FONT_SIZES['legend'],
                      framealpha=0.9, loc='upper left')
    add_panel_label(axes[0, 2], 'C')

    # ── Panel D: ROI2 ipsi vs acute behavior ──────────────────────────────────────
    rho, p = spearman_stats(ac['so_power_roi2_ipsi'].values * 1e5, y_ac)
    scatter_with_regression(
        axes[1, 0], ac['so_power_roi2_ipsi'].values * 1e5, y_ac, colors_all,
        xlabel=SO_IPSI_ROI2, ylabel=BEH_ACUTE,
        stats_text=format_stats(rho=rho, p=p), loc='upper right',
    )
    add_panel_label(axes[1, 0], 'D')

    # ── Panel E: ROI1 contra vs acute behavior (null) ─────────────────────────────
    rho, p = spearman_stats(ac['so_power_roi1_contra'].values * 1e5, y_ac)
    scatter_with_regression(
        axes[1, 1], ac['so_power_roi1_contra'].values * 1e5, y_ac, colors_all,
        xlabel=SO_CONTRA_ROI1, ylabel=BEH_ACUTE,
        stats_text=format_stats(rho=rho, p=p), loc='lower right',
    )
    add_panel_label(axes[1, 1], 'E')

    # ── Panel F: ROI2 contra vs acute behavior (null) ─────────────────────────────
    rho, p = spearman_stats(ac['so_power_roi2_contra'].values * 1e5, y_ac)
    scatter_with_regression(
        axes[1, 2], ac['so_power_roi2_contra'].values * 1e5, y_ac, colors_all,
        xlabel=SO_CONTRA_ROI2, ylabel=BEH_ACUTE,
        stats_text=format_stats(rho=rho, p=p), loc='upper right',
    )
    add_panel_label(axes[1, 2], 'F')

    plt.tight_layout(pad=0.5, h_pad=2.0, w_pad=1.5)
    save_figure(fig, 'Figure2_AcuteSO', OUTPUT_DIR)
    plt.close()
    print('Figure 2 done.')


if __name__ == '__main__':
    main()
