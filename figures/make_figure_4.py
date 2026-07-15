"""make_figure_4.py -- Figure 4: SO Recovery and Dissociation (4-panel, 1x4)

A: SO power recovery trajectory -- 4 lines (ROI1/ROI2 x ipsi/contra) % baseline, n=23
B: Bilateral SO covariation trajectory -- ROI1 and ROI2 Pearson r, n=23
C: STI fractional SO recovery -- ROI1 ipsi % baseline at week 1 (strip+box, null)
D: STI week 1 behavior -- forelimb asymmetry score (strip+box, significant)
"""

import sys
from pathlib import Path

import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import numpy as np
from scipy.stats import mannwhitneyu, pearsonr
from scipy import stats as scipy_stats

FIGURES_DIR = Path(__file__).parent
sys.path.insert(0, str(FIGURES_DIR))

from figure_style import (
    COLORS, FONT_SIZES, LINE_WIDTHS,
    set_publication_style, style_ax, add_panel_label,
    strip_box_panel, save_figure, load_data,
)

DATA_FILE = FIGURES_DIR.parent / 'WFCI AI Summary.xlsx'
OUTPUT_DIR = FIGURES_DIR / 'final'


def pearson_r(df, col1, col2):
    x, y = df[col1].values.astype(float), df[col2].values.astype(float)
    mask = np.isfinite(x) & np.isfinite(y)
    r, p = pearsonr(x[mask], y[mask])
    return r, p


def main():
    set_publication_style()
    data = load_data(DATA_FILE)
    bl, ac, wk, wk_clean = data['bl'], data['ac'], data['wk'], data['wk_clean']
    sti_pos, sti_neg = data['sti_pos'], data['sti_neg']
    sti_pos_clean = data['sti_pos_clean']
    sti_neg_clean = data['sti_neg_clean']

    # Consistent n=23 subset at all timepoints
    idx23 = wk_clean.index
    bl23  = bl.loc[idx23]
    ac23  = ac.loc[idx23]

    # Shared x-axis spec for both trajectory panels
    tp_x      = [0, 1, 2]
    tp_labels = ['Baseline\n(n=23)', 'Acute\n(24h)', 'Week 1']

    fig, axes = plt.subplots(1, 4, figsize=(18, 4.8))

    # ── Panel A: SO power recovery trajectory ────────────────────────────────
    ax = axes[0]

    line_specs_so = [
        ('so_power_roi1_ipsi',   COLORS['roi1'], '-',  'Lesion Ipsilateral'),
        ('so_power_roi1_contra', COLORS['roi1'], '--', 'Lesion Contralateral'),
        ('so_power_roi2_ipsi',   COLORS['roi2'], '-',  'Peri Ipsilateral'),
        ('so_power_roi2_contra', COLORS['roi2'], '--', 'Peri Contralateral'),
    ]

    for col, color, ls, label in line_specs_so:
        bl_vals = bl23[col].values
        ac_pct  = ac23[col].values / bl_vals * 100
        wk_pct  = wk_clean[col].values / bl_vals * 100

        means = [100.0, ac_pct.mean(), wk_pct.mean()]
        sems  = [0.0,
                 ac_pct.std(ddof=1) / np.sqrt(len(ac_pct)),
                 wk_pct.std(ddof=1) / np.sqrt(len(wk_pct))]

        ax.errorbar(tp_x, means, yerr=sems,
                    color=color, linestyle=ls,
                    lw=LINE_WIDTHS['errorbar'],
                    marker='o', markersize=3,
                    capsize=2, zorder=3, label=label)

    ax.axhline(100, color='#888888', lw=0.8, ls=':', alpha=0.5, zorder=1)
    ax.set_xticks(tp_x)
    ax.set_xticklabels(tp_labels, fontsize=FONT_SIZES['tick'])
    ax.set_ylabel('SO Power (% Baseline)', fontsize=FONT_SIZES['axis_label'])

    legend_handles_so = [
        Line2D([0], [0], color=COLORS['roi1'], lw=LINE_WIDTHS['errorbar'],
               ls='-',  marker='o', markersize=3, label='Lesion Ipsilateral'),
        Line2D([0], [0], color=COLORS['roi1'], lw=LINE_WIDTHS['errorbar'],
               ls='--', marker='o', markersize=3, label='Lesion Contralateral'),
        Line2D([0], [0], color=COLORS['roi2'], lw=LINE_WIDTHS['errorbar'],
               ls='-',  marker='o', markersize=3, label='Peri Ipsilateral'),
        Line2D([0], [0], color=COLORS['roi2'], lw=LINE_WIDTHS['errorbar'],
               ls='--', marker='o', markersize=3, label='Peri Contralateral'),
    ]
    ax.legend(handles=legend_handles_so, fontsize=4,
              loc='upper right', framealpha=0.9,
              handlelength=1.5, handleheight=0.6,
              handletextpad=0.3, borderpad=0.3, labelspacing=0.2)
    style_ax(ax)
    add_panel_label(ax, 'A')

    # ── Panel B: Bilateral SO covariation trajectory ──────────────────────────
    ax = axes[1]

    # Panel 4B uses n=25 (all animals, including Grubbs outliers DBSI_M02/DBSI_M04)
    r_bl_r1, p_bl_r1 = pearson_r(bl,  'so_power_roi1_ipsi', 'so_power_roi1_contra')
    r_ac_r1, p_ac_r1 = pearson_r(ac,  'so_power_roi1_ipsi', 'so_power_roi1_contra')
    r_wk_r1, p_wk_r1 = pearson_r(wk,  'so_power_roi1_ipsi', 'so_power_roi1_contra')

    r_bl_r2, p_bl_r2 = pearson_r(bl,  'so_power_roi2_ipsi', 'so_power_roi2_contra')
    r_ac_r2, p_ac_r2 = pearson_r(ac,  'so_power_roi2_ipsi', 'so_power_roi2_contra')
    r_wk_r2, p_wk_r2 = pearson_r(wk,  'so_power_roi2_ipsi', 'so_power_roi2_contra')

    cov_specs = [
        ([r_bl_r1, r_ac_r1, r_wk_r1], COLORS['roi1'], 'o', 'Lesion'),
        ([r_bl_r2, r_ac_r2, r_wk_r2], COLORS['roi2'], 's', 'Perilesional'),
    ]

    for r_vals, color, marker, label in cov_specs:
        ax.plot(tp_x, r_vals,
                color=color, lw=LINE_WIDTHS['errorbar'],
                marker=marker, markersize=4, zorder=3, label=label)

    ax.axhline(0, color='#888888', lw=0.8, ls='--', alpha=0.4, zorder=1)
    ax.set_xticks(tp_x)
    ax.set_xticklabels(tp_labels, fontsize=FONT_SIZES['tick'])
    ax.set_ylabel('Bilateral SO Covariation (Pearson r)',
                  fontsize=FONT_SIZES['axis_label'])
    ax.set_ylim(-0.05, 1.15)

    legend_handles_cov = [
        Line2D([0], [0], color=COLORS['roi1'], lw=LINE_WIDTHS['errorbar'],
               marker='o', markersize=4, label='Lesion'),
        Line2D([0], [0], color=COLORS['roi2'], lw=LINE_WIDTHS['errorbar'],
               marker='s', markersize=4, label='Perilesional'),
    ]
    ax.legend(handles=legend_handles_cov, fontsize=5,
              loc='lower right', framealpha=0.9,
              handlelength=2.0, handleheight=0.8,
              handletextpad=0.4, borderpad=0.4, labelspacing=0.3)
    style_ax(ax)
    add_panel_label(ax, 'B')

    # ── Panel C: STI fractional SO recovery (null) ────────────────────────────
    idx_pos = wk_clean.index[sti_pos_clean]
    idx_neg = wk_clean.index[sti_neg_clean]

    wk_so_pos = (wk_clean.loc[idx_pos, 'so_power_roi1_ipsi'].values /
                 bl.loc[idx_pos, 'so_power_roi1_ipsi'].values * 100)
    wk_so_neg = (wk_clean.loc[idx_neg, 'so_power_roi1_ipsi'].values /
                 bl.loc[idx_neg, 'so_power_roi1_ipsi'].values * 100)

    _, p_so = mannwhitneyu(wk_so_pos, wk_so_neg, alternative='two-sided')

    strip_box_panel(
        axes[2],
        data_pos=wk_so_pos,
        data_neg=wk_so_neg,
        ylabel='Lesion Ipsilateral SO\n(% Baseline, Week 1)',
        panel_label='C',
        p_val=p_so,
    )

    # ── Panel D: STI week 1 behavior (significant) ────────────────────────────
    wk_beh_pos = wk.loc[sti_pos, 'behavior_score'].values
    wk_beh_neg = wk.loc[sti_neg, 'behavior_score'].values
    _, p_beh = mannwhitneyu(wk_beh_pos, wk_beh_neg, alternative='two-sided')

    strip_box_panel(
        axes[3],
        data_pos=wk_beh_pos,
        data_neg=wk_beh_neg,
        ylabel='Week 1 Forelimb Asymmetry Score (%)',
        panel_label='D',
        p_val=p_beh,
        zero_line=True,
    )

    plt.tight_layout(pad=0.5)
    save_figure(fig, 'Figure4_RecoveryDissociation', OUTPUT_DIR)
    plt.close(fig)
    print('Figure 4 complete.')


if __name__ == '__main__':
    main()
