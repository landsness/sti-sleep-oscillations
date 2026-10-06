"""make_suppl_table_sensitivity.py -- Supplementary sensitivity table.

Recomputes every week-1-dependent statistic reported in the manuscript both
WITH the two protocol-deviation animals excluded (n = 23, the primary analysis
for week-1 SO measures) and WITH all animals retained (n = 25), to show that
no reported conclusion depends on the exclusion.

Two animals (DBSI_M02, DBSI_M04) received a structural MRI under ~1 h deep
isoflurane between the acute and week-1 imaging sessions; their week-1 SO
measurements are therefore excluded from week-1 SO analyses (see Methods).
Figure 5 (baseline LI vs week-1 behavior) uses all 25 animals as its primary
analysis because it does not depend on their week-1 SO; the n = 23 column is
its sensitivity check.

Reuses figure_style.load_data so the inputs are identical to the figures.
Outputs: results/tables/SupplTable1_Sensitivity.csv and a rendered table
(results/figures/{png,pdf,svg}/SupplTable1_Sensitivity.*). Analysis only -- no figures or
manuscript files are modified.
"""

import sys
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.stats import (wilcoxon, pearsonr, mannwhitneyu, spearmanr,
                         rankdata, linregress, norm)

FIGURES_DIR = Path(__file__).resolve().parent
sys.path[:0] = [str(FIGURES_DIR), str(FIGURES_DIR.parent), str(FIGURES_DIR.parent / 'analysis')]
import config  # noqa: E402

from figure_style import set_publication_style, save_figure, load_data

DATA_FILE = config.XLSX
OUT_DIR = config.FIG_OUT

OUTLIERS = ['DBSI_M02', 'DBSI_M04']
CHANNELS = {
    'lesion ipsi':  'so_power_roi1_ipsi',
    'lesion contra': 'so_power_roi1_contra',
    'peri ipsi':    'so_power_roi2_ipsi',
    'peri contra':  'so_power_roi2_contra',
}


# ---------------------------------------------------------------------------
# statistics helpers (identical methods to the figure scripts)
# ---------------------------------------------------------------------------
def pstars(p):
    return '***' if p < .001 else '**' if p < .01 else '*' if p < .05 else 'ns'


def pfmt(p):
    return ' < 0.001' if p < 0.001 else f' = {p:.3f}'   # includes the relation sign


def partial_spear(x, y, z):
    """Partial Spearman of x~y controlling for z (rank residuals) -- matches make_figure_5."""
    x, y, z = map(lambda a: np.asarray(a, float), (x, y, z))
    m = ~(np.isnan(x) | np.isnan(y) | np.isnan(z))
    rx, ry, rz = rankdata(x[m]), rankdata(y[m]), rankdata(z[m])

    def resid(a, b):
        sl, ic, *_ = linregress(b, a)
        return a - (sl * b + ic)
    return pearsonr(resid(rx, rz), resid(ry, rz))


def steiger1980(r12, r34, r13, r14, r23, r24, n):
    """Steiger (1980) test comparing dependent, non-overlapping correlations
    r12 (acute ipsi-contra) vs r34 (week-1 ipsi-contra) on the same animals."""
    z12, z34 = np.arctanh(r12), np.arctanh(r34)
    rbar = (r12 + r34) / 2
    num = (0.5 * r12 * r34 * (r13**2 + r14**2 + r23**2 + r24**2)
           + r13 * r24 + r14 * r23
           - (r12 * r13 * r14 + r12 * r23 * r24 + r13 * r23 * r34 + r14 * r24 * r34))
    c = num / ((1 - rbar**2)**2)
    Z = (z12 - z34) * np.sqrt((n - 3) / (2 - 2 * c))
    return Z, 2 * norm.sf(abs(Z))


# ---------------------------------------------------------------------------
# per-metric computations, each returning a formatted result string for a cohort
# ---------------------------------------------------------------------------
def traj_acute_to_wk1(bl, ac, wk, idx, col):
    b = bl.loc[idx, col].values
    ap = ac.loc[idx, col].values / b * 100
    wp = wk.loc[idx, col].values / b * 100
    _, p = wilcoxon(ap, wp)
    return f'{ap.mean():.0f}%→{wp.mean():.0f}% BL, p{pfmt(p)} {pstars(p)}', p


def traj_wk1_vs_bl(bl, wk, idx, col):
    b = bl.loc[idx, col].values
    wp = wk.loc[idx, col].values / b * 100
    _, p = wilcoxon(wp - 100)
    return f'{wp.mean():.0f}% BL, p{pfmt(p)} {pstars(p)}', p


def paired_channels(bl, wk, idx, col_a, col_b):
    a = wk.loc[idx, col_a].values / bl.loc[idx, col_a].values * 100
    b = wk.loc[idx, col_b].values / bl.loc[idx, col_b].values * 100
    _, p = wilcoxon(a, b)
    return f'{a.mean():.0f}% vs {b.mean():.0f}% BL, p{pfmt(p)} {pstars(p)}', p


def covariation_steiger(ac, wk, idx, roi):
    R = lambda u, v: pearsonr(u, v)[0]
    ai, ac_ = ac.loc[idx, f'so_power_{roi}_ipsi'], ac.loc[idx, f'so_power_{roi}_contra']
    wi, wc = wk.loc[idx, f'so_power_{roi}_ipsi'], wk.loc[idx, f'so_power_{roi}_contra']
    r12, r34 = R(ai, ac_), R(wi, wc)
    Z, p = steiger1980(r12, r34, R(ai, wi), R(ai, wc), R(ac_, wi), R(ac_, wc), len(idx))
    return f'wk1 r={r34:.3f}; z={Z:.2f}, p{pfmt(p)} {pstars(p)}'.replace('-', '\u2212'), p


def sti_recovery(bl, wk, idx, col):
    sti = bl.loc[idx, 'secondary_thalamic_injury']
    pos = idx[(sti == 1).values]
    neg = idx[(sti == 0).values]
    vp = wk.loc[pos, col].values / bl.loc[pos, col].values * 100
    vn = wk.loc[neg, col].values / bl.loc[neg, col].values * 100
    _, p = mannwhitneyu(vp, vn, alternative='two-sided')
    return f'{np.median(vp):.0f}% vs {np.median(vn):.0f}% BL (median), p{pfmt(p)} {pstars(p)}', p


def so_vs_behavior(bl, wk, idx, col):
    # week-1 SO expressed as % of individual baseline (matches Supplementary Figure 6 / former Table 1)
    pct = wk.loc[idx, col].values / bl.loc[idx, col].values * 100
    rho, p = spearmanr(pct, wk.loc[idx, 'behavior_score'].values)
    return f'ρ={rho:+.3f}, p{pfmt(p)} {pstars(p)}'.replace('-', '\u2212'), p


def baselineLI_vs_behavior(bl, wk, idx, partial=False):
    x = bl.loc[idx, 'LI_roi1'].values
    y = wk.loc[idx, 'behavior_score'].values
    if partial:
        rho, p = partial_spear(x, y, bl.loc[idx, 'stroke_size'].values)
    else:
        rho, p = spearmanr(x, y)
    return f'ρ={rho:+.3f}, p{pfmt(p)} {pstars(p)}'.replace('-', '\u2212'), p


# ---------------------------------------------------------------------------
# build table
# ---------------------------------------------------------------------------
def build_rows(data):
    bl, ac, wk = data['bl'], data['ac'], data['wk']
    i25 = wk.index
    i23 = wk.drop([o for o in OUTLIERS if o in wk.index]).index

    rows = []

    def add(section, analysis, fn23, fn25, primary):
        s23, p23 = fn23()
        s25, p25 = fn25()
        consistent = 'yes' if (p23 < .05) == (p25 < .05) else 'NO'
        rows.append(dict(Section=section, Analysis=analysis,
                         **{'Primary n': primary,
                            'n = 23 (excluded)': s23,
                            'n = 25 (all animals)': s25,
                            'Consistent': consistent}))

    # Fig 4A -- SO recovery trajectory
    for name, col in CHANNELS.items():
        add('Fig 4A', f'{name}: acute→wk1',
            lambda c=col: traj_acute_to_wk1(bl, ac, wk, i23, c),
            lambda c=col: traj_acute_to_wk1(bl, ac, wk, i25, c), 23)
    for name in ('lesion ipsi', 'peri ipsi'):
        col = CHANNELS[name]
        add('Fig 4A', f'{name}: wk1 vs baseline',
            lambda c=col: traj_wk1_vs_bl(bl, wk, i23, c),
            lambda c=col: traj_wk1_vs_bl(bl, wk, i25, c), 23)
    add('Fig 4A', 'lesion vs peri ipsi recovery',
        lambda: paired_channels(bl, wk, i23, CHANNELS['lesion ipsi'], CHANNELS['peri ipsi']),
        lambda: paired_channels(bl, wk, i25, CHANNELS['lesion ipsi'], CHANNELS['peri ipsi']), 23)
    add('Fig 4A', 'lesion wk1 ipsi vs contra',
        lambda: paired_channels(bl, wk, i23, CHANNELS['lesion ipsi'], CHANNELS['lesion contra']),
        lambda: paired_channels(bl, wk, i25, CHANNELS['lesion ipsi'], CHANNELS['lesion contra']), 23)

    # Fig 4B -- bilateral covariation (Steiger acute->wk1)
    for name, roi in (('lesion', 'roi1'), ('peri', 'roi2')):
        add('Fig 4B', f'{name} covariation acute→wk1',
            lambda r=roi: covariation_steiger(ac, wk, i23, r),
            lambda r=roi: covariation_steiger(ac, wk, i25, r), 23)

    # Fig 4C -- STI SO recovery
    add('Fig 4C', 'STI+ vs STI− lesion-ipsi recovery',
        lambda: sti_recovery(bl, wk, i23, CHANNELS['lesion ipsi']),
        lambda: sti_recovery(bl, wk, i25, CHANNELS['lesion ipsi']), 23)

    # Suppl Fig 6 (formerly Table 1) -- week-1 SO vs behavior
    for name, col in CHANNELS.items():
        add('Suppl Fig 6', f'wk1 {name} SO vs behavior',
            lambda c=col: so_vs_behavior(bl, wk, i23, c),
            lambda c=col: so_vs_behavior(bl, wk, i25, c), 23)

    # Fig 5 -- baseline LI vs week-1 behavior (primary = n=25)
    add('Fig 5', 'baseline lesion LI → wk1 behavior',
        lambda: baselineLI_vs_behavior(bl, wk, i23),
        lambda: baselineLI_vs_behavior(bl, wk, i25), 25)
    add('Fig 5', '...partial (control infarct volume)',
        lambda: baselineLI_vs_behavior(bl, wk, i23, partial=True),
        lambda: baselineLI_vs_behavior(bl, wk, i25, partial=True), 25)

    return pd.DataFrame(rows)


def render_table(df, stem):
    set_publication_style()
    cols = ['Section', 'Analysis', 'Primary n', 'n = 23 (excluded)',
            'n = 25 (all animals)', 'Consistent']
    d = df[cols]
    nrow = len(d)
    fig, ax = plt.subplots(figsize=(13, 0.34 * nrow + 1.2))
    ax.axis('off')
    tbl = ax.table(cellText=d.values, colLabels=cols, cellLoc='left', loc='center')
    tbl.auto_set_font_size(False)
    tbl.set_fontsize(7.5)
    tbl.scale(1, 1.25)
    widths = [0.07, 0.26, 0.07, 0.24, 0.24, 0.09]
    for (r, c), cell in tbl.get_celld().items():
        cell.set_edgecolor('#cccccc')
        cell.set_width(widths[c])
        if r == 0:
            cell.set_text_props(weight='bold', color='white')
            cell.set_facecolor('#4a5568')
        else:
            cell.set_facecolor('#ffffff' if r % 2 else '#f2f4f7')
            if cols[c] == 'Consistent' and d.iloc[r - 1]['Consistent'] == 'NO':
                cell.set_text_props(weight='bold', color='#c0392b')
    ax.set_title('Supplementary Table 1. Sensitivity of week-1 findings to exclusion '
                 'of two protocol-deviation animals (DBSI_M02, DBSI_M04)',
                 fontsize=8, loc='left', pad=10)
    save_figure(fig, stem, OUT_DIR)
    plt.close(fig)


def main():
    data = load_data(DATA_FILE)
    df = build_rows(data)

    config.TABLES.mkdir(parents=True, exist_ok=True)
    csv_path = config.TABLES / 'SupplTable1_Sensitivity.csv'
    df.to_csv(csv_path, index=False)

    with pd.option_context('display.max_colwidth', None, 'display.width', 200):
        print(df.to_string(index=False))
    n_incons = (df['Consistent'] == 'NO').sum()
    print(f'\nRows: {len(df)} | inconsistent conclusions: {n_incons}')
    print(f'CSV: {csv_path}')

    render_table(df, 'SupplTable1_Sensitivity')
    print('Rendered SupplTable1_Sensitivity.{png,pdf,svg}')


if __name__ == '__main__':
    main()
