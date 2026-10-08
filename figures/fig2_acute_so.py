"""fig2_acute_so.py -- Figure 2: Acute SO suppression indexes injury severity and concurrent deficit.

Rev. 2026-10-04: every panel is drawn natively at the final print size (180 mm wide,
Frontiers two-column) with one font scale (FIG2_FONTS); panels A-C are no longer
pasted in as images from figures/panels/ (their draw() functions are called here).

Layout (4 rows)
  Row 1  A: bilateral group map of acute SO power, % baseline   (panels/fig2a_lesion_topo.draw)
         B: representative lesional traces, Ms8, uncorrected vs corrected dF/F (panels/fig2b_lesional_traces.draw)
  Row 2  C: ROI power spectra, baseline vs acute, n = 25        (panels/fig2c_roi_spectra.draw)
  Row 3  D: acute SO (% baseline), 4 region x hemisphere channels; Wilcoxon brackets
         E: acute SO (% baseline) by STI group; Mann-Whitney brackets
  Row 4  F-I: acute SO power vs acute forelimb asymmetry (Spearman; BH-FDR q)
         F lesional ipsi, G perilesional ipsi, H lesional contra, I perilesional contra
"""

import sys
from pathlib import Path

import matplotlib.patches as mpatches
from matplotlib.ticker import MaxNLocator
import matplotlib.pyplot as plt
import numpy as np
from scipy import stats as scipy_stats
from scipy.stats import mannwhitneyu, wilcoxon

FIGURES_DIR = Path(__file__).resolve().parent
sys.path[:0] = [str(FIGURES_DIR), str(FIGURES_DIR / 'panels'), str(FIGURES_DIR.parent), str(FIGURES_DIR.parent / 'analysis')]
import config  # noqa: E402

from figure_style import (
    COLORS, FONT_SIZES, LINE_WIDTHS, MARKER_SIZES,
    add_panel_label, format_stats, format_pval, load_data,
    save_figure, scatter_with_regression, set_publication_style,
    style_ax, q_value,
)

DATA_FILE = config.XLSX
OUTPUT_DIR = config.FIG_OUT

# Channel order: ipsilateral left, contralateral right
COLS = ['so_power_roi1_ipsi', 'so_power_roi2_ipsi',
        'so_power_roi1_contra', 'so_power_roi2_contra']
CH_LABELS = ['Lesion', 'Peri-\nlesional', 'Lesion', 'Peri-\nlesional']
X_POS = np.array([0, 1.15, 2.75, 3.9])          # lesion/peri pairs, ipsi then contra
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


HEMI_STYLE = 'bracket'          # restored 2026-10-08 per Reviewer 1 R2: label Ipsilesional/Contralesional on D and E


def hemi_groups(ax):
    """Group the lesion/perilesional ticks by hemisphere: a thin bracket under each pair
    with the hemisphere name directly beneath it (drawn in axes-fraction y, below the ticks)."""
    if HEMI_STYLE == 'none':
        return
    import matplotlib.transforms as mtrans
    tr = mtrans.blended_transform_factory(ax.transData, ax.transAxes)
    fig = ax.figure
    # place the bracket just under the two-line tick labels (offset in points -> axes fraction)
    h_in = ax.get_position().height * fig.get_figheight()
    off_pt = 2 * 1.2 * FONT_SIZES['tick'] + 6
    yb = -off_pt / 72 / h_in
    tick = 2.5 / 72 / h_in
    for (a, b), txt in (((X_POS[0], X_POS[1]), 'Ipsilesional'), ((X_POS[2], X_POS[3]), 'Contralesional')):
        pad = 0.25
        ax.plot([a - pad, a - pad, b + pad, b + pad], [yb + tick, yb, yb, yb + tick], transform=tr,
                color='0.35', lw=0.6, clip_on=False)
        ax.text((a + b) / 2, yb - 1.5 / 72 / h_in, txt, transform=tr, ha='center', va='top',
                fontsize=FONT_SIZES['tick'], color='0.2', clip_on=False)


FIG_W = 7.09                     # in; Frontiers two-column (180 mm)
FIG_H = 7.8
FIG2_FONTS = dict(panel_label=10, title=8, axis_label=8, tick=8, stats=8, legend=8, sig_star=8)
FIG2_LINES = dict(regression=0.8, identity=0.5, axis=0.6, bracket=0.7, errorbar=0.8, median=1.2, box_whisker=0.6)
FIG2_MARKERS = dict(scatter=16, errorbar=4, strip=12)
PANEL_FS = dict(label=8, tick=8, legend=7.5, title=8, note=7)      # A-C internal text


def _label(ax, letter, fig, dx=-0.30, dy=0.08):
    """Panel letter at a fixed offset (inches) from the axes' top-left corner."""
    import matplotlib.transforms as mtrans
    off = mtrans.ScaledTranslation(dx, dy, fig.dpi_scale_trans)
    ax.text(0, 1, letter, transform=ax.transAxes + off, fontsize=FIG2_FONTS['panel_label'],
            fontweight='bold', va='bottom', ha='left')


def main():
    import figure_style as fsty
    saved = (dict(fsty.FONT_SIZES), dict(fsty.LINE_WIDTHS), dict(fsty.MARKER_SIZES))
    fsty.FONT_SIZES.update(FIG2_FONTS); fsty.LINE_WIDTHS.update(FIG2_LINES); fsty.MARKER_SIZES.update(FIG2_MARKERS)
    try:
        _build()
    finally:                                   # never leak Figure-2 sizes into other figures (make_all runs in one process)
        for d, v in zip((fsty.FONT_SIZES, fsty.LINE_WIDTHS, fsty.MARKER_SIZES), saved):
            d.clear(); d.update(v)


def _build():
    import fig2a_lesion_topo as p2a
    import fig2b_lesional_traces as p2b
    import fig2c_roi_spectra as p2c
    set_publication_style()
    plt.rcParams.update({'axes.linewidth': 0.6, 'xtick.major.width': 0.6, 'ytick.major.width': 0.6,
                         'xtick.major.size': 2.5, 'ytick.major.size': 2.5,
                         'xtick.major.pad': 1.5, 'ytick.major.pad': 1.5, 'axes.labelpad': 2})
    data = load_data(DATA_FILE)
    bl, ac = data['bl'], data['ac']
    sti_pos, sti_neg = data['sti_pos'], data['sti_neg']
    sti_pos_idx = bl.index[sti_pos]
    sti_neg_idx = bl.index[sti_neg]
    colors_all = [COLORS['sti_pos'] if sti_pos.get(mid, False) else COLORS['sti_neg'] for mid in ac.index]
    y_ac = ac['behavior_score'].values
    pct = {c: (ac[c] / bl[c] * 100).values for c in COLS}
    patch_pos = mpatches.Patch(color=COLORS['sti_pos'], label='STI+')
    patch_neg = mpatches.Patch(color=COLORS['sti_neg'], label='STI−')

    # plt.figure (not plt.subplots): this figure is already at print size, so
    # export_frontiers.py's figsize patch does not apply to it.
    fig = plt.figure(figsize=(FIG_W, FIG_H))
    outer = fig.add_gridspec(4, 1, height_ratios=[1.75, 1.35, 1.95, 1.62], hspace=0.55,
                             left=0.075, right=0.985, top=0.975, bottom=0.075)
    r1 = outer[0].subgridspec(1, 2, width_ratios=[1.0, 1.40], wspace=0.42)
    r2 = outer[1].subgridspec(1, 3, wspace=0.10)
    r3 = outer[2].subgridspec(1, 2, wspace=0.22)
    r4 = outer[3].subgridspec(2, 4, height_ratios=[0.02, 1.0], hspace=0.0, wspace=0.30)   # thin top strip = spacing below D/E tick labels

    # ── A: bilateral acute map ────────────────────────────────────────────────
    axA = fig.add_subplot(r1[0, 0])
    groups, ns, geom, outline = p2a.load_maps()
    im = p2a.draw(axA, groups['acute'], geom, outline, title=None, legend=True,
                  fs=dict(title=8, hemi=7.5, note=7, scale=7, legend=7), ring_lw=1.0)
    cb = fig.colorbar(im, ax=axA, fraction=0.045, pad=0.02, aspect=18)
    cb.set_label('SO power (% baseline)', fontsize=FIG2_FONTS['axis_label'], labelpad=2)
    cb.ax.tick_params(labelsize=FIG2_FONTS['tick'], length=2, pad=1)
    cb.outline.set_linewidth(0.5)
    _label(axA, 'A', fig, dx=-0.12)

    # ── B: representative lesional traces ─────────────────────────────────────
    rB = r1[0, 1].subgridspec(2, 1, hspace=0.12)
    axB = [fig.add_subplot(rB[0])]
    axB.append(fig.add_subplot(rB[1], sharex=axB[0], sharey=axB[0]))
    p2b.draw(axB, fs=PANEL_FS, lw=(0.6, 0.8))
    _label(axB[0], 'B', fig)

    # ── C: ROI power spectra ──────────────────────────────────────────────────
    axC = [fig.add_subplot(r2[0, 0])]
    axC += [fig.add_subplot(r2[0, k], sharey=axC[0]) for k in (1, 2)]
    p2c.draw(axC, fs=PANEL_FS, lw=1.1)
    for a in axC[1:]:
        a.tick_params(labelleft=False)
    _label(axC[0], 'C', fig)

    # ── D: all animals, 4 channels ────────────────────────────────────────────
    ax = fig.add_subplot(r3[0, 0])
    for xi, col, color in zip(X_POS, COLS, CH_COLORS):
        strip_box(ax, xi, pct[col], color)
    ax.axhline(100, color='#888888', lw=0.7, ls='--', alpha=0.6)
    _, p_r1_hemi = wilcoxon(pct['so_power_roi1_ipsi'], pct['so_power_roi1_contra'])
    _, p_r2_hemi = wilcoxon(pct['so_power_roi2_ipsi'], pct['so_power_roi2_contra'])
    _, p_gradient = wilcoxon(pct['so_power_roi1_ipsi'], pct['so_power_roi2_ipsi'])
    ymax = max(pct[c].max() for c in COLS)
    sig_bracket(ax, X_POS[0], X_POS[1], ymax + 8, p_gradient, dy=3)
    sig_bracket(ax, X_POS[0], X_POS[2], ymax + 30, p_r1_hemi, dy=3)
    sig_bracket(ax, X_POS[1], X_POS[3], ymax + 52, p_r2_hemi, dy=3)
    ax.set_xticks(X_POS); ax.set_xticklabels(CH_LABELS)
    hemi_groups(ax)
    ax.set_ylabel('SO power (% baseline)')
    ax.set_ylim(-5, ymax + 72)
    style_ax(ax)
    _label(ax, 'D', fig)

    # ── E: STI+ vs STI− ───────────────────────────────────────────────────────
    ax = fig.add_subplot(r3[0, 1])
    w = 0.30
    prev_y = -np.inf
    for xi, col in zip(X_POS, COLS):
        vpos = (ac.loc[sti_pos_idx, col] / bl.loc[sti_pos_idx, col] * 100).values
        vneg = (ac.loc[sti_neg_idx, col] / bl.loc[sti_neg_idx, col] * 100).values
        strip_box(ax, xi - w * 0.65, vpos, COLORS['sti_pos'], width=w)
        strip_box(ax, xi + w * 0.65, vneg, COLORS['sti_neg'], width=w)
        _, p_mw = mannwhitneyu(vpos, vneg, alternative='two-sided')
        y_top = max(vpos.max(), vneg.max()) + 6
        if y_top < prev_y + 16:          # neighbouring p-labels would collide -> raise this one
            y_top = prev_y + 16
        sig_bracket(ax, xi - w * 0.65, xi + w * 0.65, y_top, p_mw, dy=3)
        prev_y = y_top
    ax.axhline(100, color='#888888', lw=0.7, ls='--', alpha=0.6)
    ax.set_xticks(X_POS); ax.set_xticklabels(CH_LABELS)
    hemi_groups(ax)
    ax.set_ylabel('SO power (% baseline)')
    ax.set_ylim(-5, prev_y + 22)
    ax.legend(handles=[patch_pos, patch_neg], loc='upper left', frameon=False,
              handlelength=1.0, handleheight=0.8, borderaxespad=0.2)
    style_ax(ax)
    _label(ax, 'E', fig)

    # ── F–I: acute SO power vs acute behaviour ────────────────────────────────
    spec = [('so_power_roi1_ipsi', 'Lesional, ipsilesional', 'lesion ipsi', 'F'),
            ('so_power_roi2_ipsi', 'Perilesional, ipsilesional', 'peri ipsi', 'G'),
            ('so_power_roi1_contra', 'Lesional, contralesional', 'lesion contra', 'H'),
            ('so_power_roi2_contra', 'Perilesional, contralesional', 'peri contra', 'I')]
    axs4 = []
    for k, (col, title, qkey, letter) in enumerate(spec):
        ax = fig.add_subplot(r4[1, k], sharey=axs4[0] if axs4 else None)
        x = ac[col].values * 1e5
        rho, p = spearman_stats(x, y_ac)
        scatter_with_regression(ax, x, y_ac, colors_all, xlabel='', ylabel='',
                                stats_text=format_stats(rho=rho, p=p, q=q_value(f'acute_so_vs_acute_behavior_{qkey}', p)),
                                loc='upper right')
        ax.set_title(title, fontsize=FIG2_FONTS['title'], pad=3)
        ax.xaxis.set_major_locator(MaxNLocator(nbins=4, integer=True))
        if k == 0:
            ax.set_ylabel('Acute forelimb\nasymmetry (%)')
        else:
            ax.tick_params(labelleft=False)
        _label(ax, letter, fig, dx=-0.30 if k == 0 else -0.16)
        axs4.append(ax)
    # headroom so the upper-right stats boxes sit in whitespace and never overlap a
    # data point (panels H, I shared the y-axis with F, so set it once on axs4[0])
    ylo, yhi = axs4[0].get_ylim()
    axs4[0].set_ylim(ylo, yhi + 0.35 * (yhi - ylo))
    # one shared x label under the row
    fig.text((axs4[0].get_position().x0 + axs4[-1].get_position().x1) / 2,
             axs4[0].get_position().y0 - 0.034,
             'Acute SO power (\u00d710$^{-5}$ (\u0394F/F)\u00b2/Hz)',
             ha='center', va='top', fontsize=FIG2_FONTS['axis_label'])

    save_figure(fig, 'Figure2_AcuteSO', OUTPUT_DIR)
    plt.close(fig)
    print('Figure 2 done.')


if __name__ == '__main__':
    main()
