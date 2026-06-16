"""make_figure_3.py — Figure 3: Bilateral SO Covariation (4-panel scatter)

Panels:
    A: Baseline ROI1  — ipsilateral vs contralateral SO power
    B: Baseline ROI2  — ipsilateral vs contralateral SO power
    C: Acute ROI1     — ipsilateral vs contralateral SO power
    D: Acute ROI2     — ipsilateral vs contralateral SO power
"""

import sys
from pathlib import Path

import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import numpy as np
from scipy.stats import pearsonr

# ---------------------------------------------------------------------------
# Path setup — allow importing figure_style from sibling directory
# ---------------------------------------------------------------------------
FIGURES_DIR = Path(__file__).parent
sys.path.insert(0, str(FIGURES_DIR))

from figure_style import (
    COLORS, FONT_SIZES, LINE_WIDTHS, MARKER_SIZES,
    set_publication_style, style_ax, add_panel_label, add_stats_box,
    scatter_with_regression, add_identity_line, format_pval, save_figure,
    load_data,
)

# ---------------------------------------------------------------------------
# Data file
# ---------------------------------------------------------------------------
DATA_FILE = FIGURES_DIR.parent / 'WFCI AI Summary.xlsx'

# ---------------------------------------------------------------------------
# Data loading and preparation
# ---------------------------------------------------------------------------

def load_data_fig3():
    return load_data(DATA_FILE)


def prepare_data(data):
    bl = data['bl']
    ac = data['ac']
    sti_pos = data['sti_pos']
    sti_neg = data['sti_neg']
    colors_all = [COLORS['sti_pos'] if sti_pos[mid] else COLORS['sti_neg'] for mid in bl.index]
    return dict(bl=bl, ac=ac, sti_pos=sti_pos, sti_neg=sti_neg, colors_all=colors_all)


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

sti_col = 'secondary_thalamic_injury'


def _point_colors(df):
    """Return list of colors based on STI group for each row in df."""
    return [
        COLORS['sti_pos'] if df.loc[idx, sti_col] == 1 else COLORS['sti_neg']
        for idx in df.index
    ]


def _compute_pearson(x, y):
    mask = np.isfinite(x) & np.isfinite(y)
    if mask.sum() < 3:
        return np.nan, np.nan
    r, p = pearsonr(x[mask], y[mask])
    return r, p


# ---------------------------------------------------------------------------
# Build figure
# ---------------------------------------------------------------------------

def make_figure(prepared):
    bl = prepared['bl']
    ac = prepared['ac']

    set_publication_style()

    fig, axes = plt.subplots(2, 2, figsize=(10, 9))

    SCALE = 1e5  # scale SO power values for display

    # Matrix layout: baseline = top row, acute = bottom row; ROI1 = left col, ROI2 = right col
    # A=BL ROI1, B=BL ROI2, C=Acute ROI1, D=Acute ROI2
    # NOTE: manuscript in-text citations must use 3A,B for baseline and 3C,D for acute
    timepoints = [
        ('Baseline',    bl, 'so_power_roi1_ipsi', 'so_power_roi1_contra', 'ROI1', 'A'),
        ('Baseline',    bl, 'so_power_roi2_ipsi', 'so_power_roi2_contra', 'ROI2', 'B'),
        ('Acute (24h)', ac, 'so_power_roi1_ipsi', 'so_power_roi1_contra', 'ROI1', 'C'),
        ('Acute (24h)', ac, 'so_power_roi2_ipsi', 'so_power_roi2_contra', 'ROI2', 'D'),
    ]

    ax_flat = [axes[0, 0], axes[0, 1], axes[1, 0], axes[1, 1]]

    # Legend patches — only added to panel A
    patch_pos = mpatches.Patch(color=COLORS['sti_pos'], label='STI+')
    patch_neg = mpatches.Patch(color=COLORS['sti_neg'], label='STI−')

    for ax, (tp_label, df, ipsi_col, contra_col, roi, panel) in zip(ax_flat, timepoints):
        x = df[ipsi_col].values.astype(float) * SCALE
        y = df[contra_col].values.astype(float) * SCALE
        colors = _point_colors(df)

        r, p = _compute_pearson(x, y)
        stats_text = f'r = {r:.3f}\n{format_pval(p)}' if not np.isnan(r) else ''

        ax.scatter(
            x, y,
            c=colors,
            s=MARKER_SIZES['scatter'],
            alpha=0.85,
            edgecolors='white',
            lw=0.5,
            zorder=3,
        )

        # Solid regression line
        mask = np.isfinite(x) & np.isfinite(y)
        if mask.sum() > 1:
            from scipy.stats import linregress
            slope, intercept, *_ = linregress(x[mask], y[mask])
            x_line = np.array([x[mask].min(), x[mask].max()])
            ax.plot(
                x_line, slope * x_line + intercept,
                color=COLORS['neutral'],
                lw=LINE_WIDTHS['regression'],
                linestyle='-',
                alpha=0.75,
                zorder=2,
            )

        style_ax(ax)

        # Identity line (dashed gray)
        add_identity_line(ax)

        # Stats annotation (upper left)
        if stats_text:
            add_stats_box(ax, stats_text, loc='upper left')

        ax.set_xlabel(f'{roi} Ipsilateral SO Power ($\\times 10^{{-5}}$ µV²)',
                      fontsize=FONT_SIZES['axis_label'])
        ax.set_ylabel(f'{roi} Contralateral SO Power ($\\times 10^{{-5}}$ µV²)',
                      fontsize=FONT_SIZES['axis_label'])
        # No panel title — handled by panel labels and figure legend

        add_panel_label(ax, panel)

        # Legend only in panel A
        if panel == 'A':
            ax.legend(
                handles=[patch_pos, patch_neg],
                fontsize=FONT_SIZES['legend'],
                framealpha=0.85,
                loc='lower right',
            )

    plt.tight_layout(pad=0.5)

    # Add column headers and row labels to make the matrix self-documenting
    col_headers = ['ROI1 (Infarct Core)', 'ROI2 (Perilesional)']
    row_labels  = ['Baseline', 'Acute (24h)']

    for col_idx, header in enumerate(col_headers):
        ax = ax_flat[col_idx]  # top row axes
        ax.annotate(header,
                    xy=(0.5, 1.0), xycoords='axes fraction',
                    xytext=(0, 28), textcoords='offset points',
                    ha='center', va='bottom',
                    fontsize=FONT_SIZES['axis_label'],
                    fontweight='bold')

    for row_idx, label in enumerate(row_labels):
        ax = ax_flat[row_idx * 2]  # left column axes
        ax.annotate(label,
                    xy=(0.0, 0.5), xycoords='axes fraction',
                    xytext=(-52, 0), textcoords='offset points',
                    ha='right', va='center',
                    fontsize=FONT_SIZES['axis_label'],
                    fontweight='bold', rotation=90)

    return fig


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    data = load_data_fig3()
    prepared = prepare_data(data)
    fig = make_figure(prepared)

    OUTPUT_DIR = FIGURES_DIR / 'final'
    save_figure(fig, 'Figure3_BilateralSOCovariation', OUTPUT_DIR)
    plt.close(fig)
    print('Figure 3 complete.')


if __name__ == '__main__':
    main()
