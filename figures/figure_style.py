"""figure_style.py — Shared plotting style library for SO/stroke manuscript."""

import os
from pathlib import Path

import matplotlib
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import numpy as np
from scipy import stats

__all__ = [
    'COLORS', 'FONT_SIZES', 'LINE_WIDTHS', 'MARKER_SIZES', 'EXPORT_SETTINGS',
    'set_publication_style', 'style_ax', 'add_panel_label', 'add_stats_box',
    'scatter_with_regression', 'add_identity_line', 'add_sig_bracket',
    'format_pval', 'format_qval', 'format_stats', 'strip_box_panel', 'coupling_bar_chart',
    'save_figure',
]

# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------

COLORS = {
    'sti_pos': '#C0392B',      # STI+ (red)
    'sti_neg': '#27AE60',      # STI- (green)
    'roi1': '#1A5276',         # ROI1 / Infarct Core (dark blue)
    'roi2': '#5DADE2',         # ROI2 / Perilesional (light blue)
    'baseline': '#2C3E50',     # Timepoint: baseline (dark navy)
    'acute': '#E67E22',        # Timepoint: acute 24h (orange)
    'week1': '#2980B9',        # Timepoint: week 1 (blue)
    'regression': '#333333',   # Regression line color
    'identity': '#AAAAAA',     # Identity line color
    'neutral': '#555555',      # Neutral gray
    'roi2_raw': '#2C3E50',     # SupplFig_S7 raw rho bar
    'roi2_partial': '#85C1E9', # SupplFig_S7 partial rho bar
}

FONT_SIZES = {
    'panel_label': 11,   # Bold A, B, C panel labels
    'title': 10,         # Axis/panel titles (removed for pub — kept for dev)
    'axis_label': 9,     # x/y axis labels
    'tick': 8,           # Tick labels
    'stats': 8,          # Stats annotation text
    'legend': 8,         # Legend text
    'sig_star': 10,      # Significance stars
}

LINE_WIDTHS = {
    'regression': 1.5,
    'identity': 1.0,
    'axis': 0.8,
    'bracket': 1.0,
    'errorbar': 1.8,
    'median': 2.0,
    'box_whisker': 1.0,
}

MARKER_SIZES = {
    'scatter': 55,    # scatter s= parameter
    'errorbar': 7,    # errorbar markersize= parameter
    'strip': 40,      # strip plot s= parameter
}

EXPORT_SETTINGS = {
    'dpi': 600,
    'svg_fonttype': 'none',
    'formats': ('svg', 'pdf', 'png'),
    'bbox_inches': 'tight',
    'pad_inches': 0.05,
}

# ---------------------------------------------------------------------------
# Global style
# ---------------------------------------------------------------------------

def set_publication_style():
    # Try preferred sans-serif fonts in order; matplotlib falls back automatically
    matplotlib.rcParams.update({
        'font.family': 'sans-serif',
        'font.sans-serif': ['Liberation Sans', 'Arimo', 'DejaVu Sans'],
        'font.size': FONT_SIZES['tick'],
        'axes.titlesize': FONT_SIZES['title'],
        'axes.labelsize': FONT_SIZES['axis_label'],
        'xtick.labelsize': FONT_SIZES['tick'],
        'ytick.labelsize': FONT_SIZES['tick'],
        'legend.fontsize': FONT_SIZES['legend'],
        'axes.spines.top': False,
        'axes.spines.right': False,
        'axes.linewidth': LINE_WIDTHS['axis'],
        'xtick.major.width': LINE_WIDTHS['axis'],
        'ytick.major.width': LINE_WIDTHS['axis'],
        'svg.fonttype': EXPORT_SETTINGS['svg_fonttype'],
        'savefig.dpi': EXPORT_SETTINGS['dpi'],
        'savefig.bbox': 'tight',
        'savefig.pad_inches': EXPORT_SETTINGS['pad_inches'],
        'figure.dpi': 100,
    })

# ---------------------------------------------------------------------------
# Axes helpers
# ---------------------------------------------------------------------------

def style_ax(ax):
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.tick_params(labelsize=FONT_SIZES['tick'])
    return ax


def add_panel_label(ax, label, x=-0.12, y=1.05):
    ax.text(
        x, y, label,
        transform=ax.transAxes,
        fontsize=FONT_SIZES['panel_label'],
        fontweight='bold',
        va='top', ha='left',
    )
    return ax


def add_stats_box(ax, text, loc='upper left'):
    # Map loc string to (x, y, ha, va) in axes coordinates
    _loc_map = {
        'upper left':  (0.05, 0.95, 'left',  'top'),
        'upper right': (0.95, 0.95, 'right', 'top'),
        'lower left':  (0.05, 0.05, 'left',  'bottom'),
        'lower right': (0.95, 0.05, 'right', 'bottom'),
    }
    x, y, ha, va = _loc_map.get(loc, _loc_map['upper left'])
    ax.text(
        x, y, text,
        transform=ax.transAxes,
        fontsize=FONT_SIZES['stats'],
        ha=ha, va=va,
        bbox=dict(boxstyle='round,pad=0.3', facecolor='white', alpha=0.85, edgecolor='none'),
    )
    return ax

# ---------------------------------------------------------------------------
# Scatter + regression
# ---------------------------------------------------------------------------

def scatter_with_regression(ax, x, y, colors, xlabel, ylabel, stats_text, loc='upper right',
                            regression=True):
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)

    # Drop NaN pairs before regression
    mask = np.isfinite(x) & np.isfinite(y)
    xc, yc = x[mask], y[mask]

    ax.scatter(
        x, y,
        c=colors,
        s=MARKER_SIZES['scatter'],
        alpha=0.85,
        edgecolors='white',
        lw=0.5,
        zorder=3,
    )

    if regression and xc.size > 1:
        slope, intercept, *_ = stats.linregress(xc, yc)
        x_line = np.array([xc.min(), xc.max()])
        ax.plot(
            x_line, slope * x_line + intercept,
            color=COLORS['regression'],
            lw=LINE_WIDTHS['regression'],
            linestyle='--',
            alpha=0.7,
            zorder=2,
        )

    style_ax(ax)
    ax.set_xlabel(xlabel, fontsize=FONT_SIZES['axis_label'])
    ax.set_ylabel(ylabel, fontsize=FONT_SIZES['axis_label'])
    add_stats_box(ax, stats_text, loc=loc)
    return ax


def add_identity_line(ax):
    # Compute limits from current view before drawing so the line spans the data range
    xlim = ax.get_xlim()
    ylim = ax.get_ylim()
    lo = min(xlim[0], ylim[0])
    hi = max(xlim[1], ylim[1])
    ax.plot(
        [lo, hi], [lo, hi],
        color=COLORS['identity'],
        lw=LINE_WIDTHS['identity'],
        linestyle='--',
        alpha=0.6,
        zorder=1,
    )
    ax.set_xlim(xlim)
    ax.set_ylim(ylim)
    return ax

# ---------------------------------------------------------------------------
# Significance bracket
# ---------------------------------------------------------------------------

def add_sig_bracket(ax, x1, x2, y, label, dy_frac=0.02):
    ylim = ax.get_ylim()
    dy = (ylim[1] - ylim[0]) * dy_frac  # tick height in data coordinates

    # Horizontal bar
    ax.plot([x1, x2], [y, y], color='black', lw=LINE_WIDTHS['bracket'])
    # Descending ticks at each end
    ax.plot([x1, x1], [y - dy, y], color='black', lw=LINE_WIDTHS['bracket'])
    ax.plot([x2, x2], [y - dy, y], color='black', lw=LINE_WIDTHS['bracket'])

    ax.text(
        (x1 + x2) / 2, y + dy * 0.5,
        label,
        ha='center', va='bottom',
        fontsize=FONT_SIZES['stats'],
    )
    return ax

# ---------------------------------------------------------------------------
# p-value / stats formatting
# ---------------------------------------------------------------------------

def format_pval(p):
    if p < 0.001:
        return 'p < 0.001'
    return f'p = {p:.3f}'


def format_qval(q):
    """Benjamini-Hochberg FDR-adjusted p (confirmatory correlation family; analysis/fdr.py)."""
    if q < 0.001:
        return 'q < 0.001'
    return f'q = {q:.3f}'


def format_stats(rho=None, r=None, p=None, n=None, partial_rho=None, partial_p=None,
                 show_n=False, q=None):
    # n is hidden by default — not shown on figures per publication style
    lines = []
    if rho is not None:
        lines.append(f'ρ = {rho:.3f}'.replace('-', '\u2212'))
    if r is not None:
        lines.append(f'r = {r:.3f}'.replace('-', '\u2212'))
    if p is not None:
        lines.append(format_pval(p))
    if q is not None:
        lines.append(format_qval(q))
    if n is not None and show_n:
        lines.append(f'n = {n}')
    if partial_rho is not None:
        lines.append(f'partial ρ = {partial_rho:.3f}'.replace('-', '\u2212'))
        if partial_p is not None:
            lines.append(format_pval(partial_p))
    return '\n'.join(lines)

# ---------------------------------------------------------------------------
# Strip + box panel
# ---------------------------------------------------------------------------

def strip_box_panel(ax, data_pos, data_neg, ylabel, panel_label, p_val=None,
                    u_val=None, sig_text='', zero_line=False):
    # u_val accepted but intentionally not displayed (kept for back-compat)
    # zero_line=True draws a dashed reference line at y=0 (e.g. forelimb asymmetry)
    data_pos = np.asarray(data_pos, dtype=float)
    data_neg = np.asarray(data_neg, dtype=float)

    np.random.seed(42)  # reproducible jitter

    bp = ax.boxplot(
        [data_pos, data_neg],
        positions=[0, 0.65],
        widths=0.28,           # narrower — publication standard
        patch_artist=True,
        notch=False,
        medianprops=dict(color='black', lw=LINE_WIDTHS['median']),
        whiskerprops=dict(lw=LINE_WIDTHS['box_whisker']),
        capprops=dict(lw=LINE_WIDTHS['box_whisker']),
        boxprops=dict(lw=LINE_WIDTHS['box_whisker']),
        flierprops=dict(marker='', linestyle='none'),
        zorder=2,
    )

    # Transparent fill so strip dots are visible through box
    for patch, color in zip(bp['boxes'], [COLORS['sti_pos'], COLORS['sti_neg']]):
        patch.set_facecolor(color)
        patch.set_alpha(0.25)

    # Jittered strip
    jitter_scale = 0.06
    for pos, data, color in zip([0, 0.65], [data_pos, data_neg], [COLORS['sti_pos'], COLORS['sti_neg']]):
        jitter = np.random.uniform(-jitter_scale, jitter_scale, size=len(data))
        ax.scatter(
            pos + jitter, data,
            color=color,
            s=MARKER_SIZES['strip'],
            alpha=0.75,
            edgecolors='white',
            lw=0.3,
            zorder=3,
        )

    ax.set_xticks([0, 0.65])
    ax.set_xticklabels(['STI+', 'STI−'], fontsize=FONT_SIZES['tick'])
    ax.set_xlim(-0.35, 1.0)   # tight margins around the two groups
    ax.set_ylabel(ylabel, fontsize=FONT_SIZES['axis_label'])
    # No panel title — description belongs in figure legend

    # Reference line at y=0 (e.g. perfect forelimb symmetry)
    if zero_line:
        ax.axhline(0, color='#666666', linewidth=0.8, linestyle='--',
                   alpha=0.6, zorder=1)

    # Significance bracket
    if p_val is not None:
        sig = p_val < 0.05
        # Position bracket above the highest data point + 8% of range
        all_data = np.concatenate([data_pos, data_neg])
        yrange = np.nanmax(all_data) - np.nanmin(all_data)
        y_bracket = np.nanmax(all_data) + 0.08 * yrange
        tick_h   = 0.025 * yrange

        # Bracket color and weight: black+bold for sig, gray+regular for ns
        bcolor = 'black' if sig else '#999999'
        blw    = 1.2 if sig else 0.9

        # Horizontal line spanning both group positions
        ax.plot([0, 0.65], [y_bracket, y_bracket], color=bcolor, lw=blw, zorder=4)
        # Vertical ticks
        ax.plot([0, 0],       [y_bracket - tick_h, y_bracket], color=bcolor, lw=blw, zorder=4)
        ax.plot([0.65, 0.65], [y_bracket - tick_h, y_bracket], color=bcolor, lw=blw, zorder=4)

        # p-value label: bold + black if significant, gray if not
        ptext = format_pval(p_val)
        fw = 'bold' if sig else 'normal'
        fc = 'black' if sig else '#999999'
        ax.text(0.325, y_bracket + 0.015 * yrange, ptext,
                ha='center', va='bottom',
                fontsize=FONT_SIZES['stats'],
                fontweight=fw, color=fc,
                zorder=5)

        # Extend y-axis to fully show bracket + label
        ax.set_ylim(top=y_bracket + 0.12 * yrange)

    style_ax(ax)
    add_panel_label(ax, panel_label)
    return ax

# ---------------------------------------------------------------------------
# Coupling bar chart
# ---------------------------------------------------------------------------

def coupling_bar_chart(ax, r_vals, p_vals, tp_labels, roi_label, panel_label, bar_color, fisher_pairs=None, fisher_color=None):
    r_vals = np.asarray(r_vals, dtype=float)
    p_vals = np.asarray(p_vals, dtype=float)
    x = np.arange(len(r_vals))

    ax.bar(x, r_vals, color=bar_color, width=0.55, zorder=2)

    # r= annotation above each bar (or below if negative)
    for xi, rv, pv in zip(x, r_vals, p_vals):
        offset = 0.02 if rv >= 0 else -0.04
        va = 'bottom' if rv >= 0 else 'top'
        ax.text(
            xi, rv + offset,
            f'r={rv:.2f}',
            ha='center', va=va,
            fontsize=FONT_SIZES['stats'],
        )
        # Significance stars above r label
        star = _sig_stars(pv)
        if star:
            star_offset = offset + (0.06 if rv >= 0 else -0.06)
            ax.text(
                xi, rv + star_offset,
                star,
                ha='center', va=va,
                fontsize=FONT_SIZES['sig_star'],
            )

    # Fisher z comparison brackets
    if fisher_pairs:
        y_top = np.nanmax(np.abs(r_vals)) + 0.15  # start brackets above bars
        fc = fisher_color if fisher_color else COLORS['neutral']
        for i, j, z_stat, fp in fisher_pairs:
            star = _sig_stars(fp)
            label = f'z={z_stat:.2f}\n{star}' if star else f'z={z_stat:.2f}\n{format_pval(fp)}'
            add_sig_bracket(ax, i, j, y_top, label)
            y_top += 0.12  # stack additional brackets

    ax.set_xticks(x)
    ax.set_xticklabels(tp_labels, fontsize=FONT_SIZES['tick'])
    ax.set_ylabel(f'Pearson r ({roi_label})', fontsize=FONT_SIZES['axis_label'])
    ax.axhline(0, color='black', lw=LINE_WIDTHS['axis'], zorder=1)

    style_ax(ax)
    add_panel_label(ax, panel_label)
    return ax


def _sig_stars(p):
    # Internal helper; not exported
    if p < 0.001:
        return '***'
    if p < 0.01:
        return '**'
    if p < 0.05:
        return '*'
    return ''

# ---------------------------------------------------------------------------
# Figure export
# ---------------------------------------------------------------------------

def save_figure(fig, output_stem, outdir, formats=('svg', 'pdf', 'png')):
    outdir = Path(outdir)
    paths = {}
    for fmt in formats:
        subdir = outdir / fmt
        subdir.mkdir(parents=True, exist_ok=True)
        fpath = subdir / f'{output_stem}.{fmt}'
        fig.savefig(
            fpath,
            format=fmt,
            dpi=EXPORT_SETTINGS['dpi'],
            bbox_inches=EXPORT_SETTINGS['bbox_inches'],
            pad_inches=EXPORT_SETTINGS['pad_inches'],
        )
        paths[fmt] = str(fpath)
        print(f'Saved: {fpath}')
    return paths

# ---------------------------------------------------------------------------
# Shared data loader -- lives in analysis/common.py; re-exported here so that
# every figure script keeps a single import line.
# ---------------------------------------------------------------------------
import sys as _sys
_ROOT = Path(__file__).resolve().parents[1]
for _p in (str(_ROOT), str(_ROOT / 'analysis')):
    if _p not in _sys.path:
        _sys.path.insert(0, _p)
from common import load_data, partial_spearman, steiger_dependent  # noqa: E402,F401
from fdr import q_value  # noqa: E402,F401


# ---------------------------------------------------------------------------
# Composite figures: raster/vector panels built by figures/panels/*.py are
# placed into Figures 1 and 2 as images.
# ---------------------------------------------------------------------------
def panel_image(ax, stem, builder=None):
    """Show results/panels/<stem>.png in `ax` (building it first with the
    panel script `builder` -- a module name in figures/panels -- if missing)."""
    import importlib
    import config
    png = config.PANELS / f'{stem}.png'
    if not png.exists() and builder is not None:
        panels_dir = str(Path(__file__).resolve().parent / 'panels')
        if panels_dir not in _sys.path:
            _sys.path.insert(0, panels_dir)
        mod = importlib.import_module(builder)
        if hasattr(mod, 'main'):
            mod.main()
    img = plt.imread(png)
    ax.imshow(img, interpolation='lanczos')
    ax.set_axis_off()
    return img.shape[1] / img.shape[0]          # aspect (w / h)


def image_row(fig, gs_row, stems_builders, width_ratios=None):
    """Fill one GridSpec row with image panels; returns their axes."""
    import config
    if width_ratios is None:
        width_ratios = []
        for stem, _ in stems_builders:
            png = config.PANELS / f'{stem}.png'
            try:
                h, w = plt.imread(png).shape[:2]
                width_ratios.append(w / h)
            except FileNotFoundError:
                width_ratios.append(1.0)
    sub = gs_row.subgridspec(1, len(stems_builders), width_ratios=width_ratios, wspace=0.04)
    axes = []
    for i, (stem, builder) in enumerate(stems_builders):
        ax = fig.add_subplot(sub[0, i])
        panel_image(ax, stem, builder)
        axes.append(ax)
    return axes
