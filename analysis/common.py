"""common.py -- shared data loading and statistical helpers.

Every figure, table and statistic in the manuscript is computed from the
DataFrames returned by `load_data()`, so all scripts use identical cohort
definitions, laterality index and %-baseline normalisation.
"""
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import pearsonr, rankdata, linregress, norm

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import config  # noqa: E402


def load_data(data_file=None):
    """Load the curated WFCI table and return per-timepoint DataFrames.

    Returns dict with keys: bl, ac, wk, wk_clean, sti_pos, sti_neg,
    sti_pos_clean, sti_neg_clean. All DataFrames are indexed by mouse_id.
      * wk_clean drops the two animals whose week-1 SO is excluded (n = 23)
      * LI_<roi> = (ipsi - contra) / (ipsi + contra)
      * <so_power col>_pct = timepoint / own baseline * 100
    """
    data_file = config.XLSX if data_file is None else data_file
    df = pd.read_excel(data_file, sheet_name=config.XLSX_SHEET)
    df.columns = df.columns.str.strip()
    for col in df.select_dtypes(include=['object', 'string']).columns:
        df[col] = df[col].astype(str).str.strip()

    bl = df[df['time_point'] == 'baseline'].set_index('mouse_id').copy()
    ac = df[df['time_point'] == '24 hours'].set_index('mouse_id').copy()
    wk = df[df['time_point'] == '1 week'].set_index('mouse_id').copy()

    wk_clean = wk.drop([o for o in config.WEEK1_SO_EXCLUDED if o in wk.index])

    sti_pos = bl['secondary_thalamic_injury'] == 1
    sti_neg = bl['secondary_thalamic_injury'] == 0
    bl_clean_sti = bl.loc[wk_clean.index, 'secondary_thalamic_injury']
    sti_pos_clean = bl_clean_sti == 1
    sti_neg_clean = bl_clean_sti == 0

    for tp in [bl, ac, wk, wk_clean]:
        for roi in ['roi1', 'roi2']:
            ipsi = tp[f'so_power_{roi}_ipsi']
            contra = tp[f'so_power_{roi}_contra']
            tp[f'LI_{roi}'] = (ipsi - contra) / (ipsi + contra)
        for col in [c for c in tp.columns if c.startswith('so_power_')]:
            if col in bl.columns:
                tp[col + '_pct'] = (tp[col] / bl.loc[tp.index, col]) * 100

    return dict(bl=bl, ac=ac, wk=wk, wk_clean=wk_clean,
                sti_pos=sti_pos, sti_neg=sti_neg,
                sti_pos_clean=sti_pos_clean, sti_neg_clean=sti_neg_clean)


def partial_spearman(x, y, z):
    """Partial Spearman rho of x~y controlling for z (rank-residual method).
    Ranks x, y, z; regresses the ranks of x and y on the ranks of z; returns
    Pearson r (and p) of the residuals."""
    x, y, z = (np.asarray(a, float) for a in (x, y, z))
    m = ~(np.isnan(x) | np.isnan(y) | np.isnan(z))
    rx, ry, rz = rankdata(x[m]), rankdata(y[m]), rankdata(z[m])

    def resid(a, b):
        fit = linregress(b, a)
        return a - (fit.slope * b + fit.intercept)
    return pearsonr(resid(rx, rz), resid(ry, rz))


def steiger_dependent(r12, r34, r13, r14, r23, r24, n):
    """Test of two dependent, non-overlapping correlations r12 vs r34 measured
    on the same n subjects (Steiger 1980). This is the implementation behind
    the manuscript's acute->week-1 covariation values; see
    results/stats for the method audit."""
    z12, z34 = np.arctanh(r12), np.arctanh(r34)
    rbar = (r12 + r34) / 2
    num = (0.5 * r12 * r34 * (r13**2 + r14**2 + r23**2 + r24**2)
           + r13 * r24 + r14 * r23
           - (r12 * r13 * r14 + r12 * r23 * r24 + r13 * r23 * r34 + r14 * r24 * r34))
    c = num / ((1 - rbar**2) ** 2)
    Z = (z12 - z34) * np.sqrt((n - 3) / (2 - 2 * c))
    return Z, 2 * norm.sf(abs(Z))


def roi_mouse_id(m):
    """ROI_Summary_V2 'Mouse' value (1, 2, ..., '2DBSI') -> analysis mouse_id ('M01', 'DBSI_M02')."""
    m = str(m)
    return f'DBSI_M0{m[0]}' if 'DBSI' in m else f'M{int(m):02d}'


# ---------------------------------------------------------------------------
# Effect sizes and 95% confidence intervals (Reviewer 2, minor point 2)
# ---------------------------------------------------------------------------
BOOT_N, BOOT_SEED = 10000, 20260929


def corr_ci(r, n, kind='spearman', n_covariates=0, level=0.95):
    """Fisher-z confidence interval for a correlation coefficient.
    Spearman rho uses the Bonett & Wright (2000) standard error
    sqrt((1 + rho^2/2) / (n - 3 - k)); Pearson r uses 1 / sqrt(n - 3 - k).
    k = number of partialled covariates."""
    df = n - 3 - n_covariates
    se = np.sqrt((1 + r ** 2 / 2) / df) if kind == 'spearman' else 1 / np.sqrt(df)
    zc = norm.ppf(0.5 + level / 2)
    z = np.arctanh(r)
    return float(np.tanh(z - zc * se)), float(np.tanh(z + zc * se))


def _rrb_ind(x, y):
    """Rank-biserial r for two independent samples: P(x > y) - P(x < y) (Kerby 2014)."""
    x, y = np.asarray(x, float), np.asarray(y, float)
    diff = x[:, None] - y[None, :]
    return float((np.sign(diff)).mean())


def _rrb_paired(d):
    """Matched-pairs rank-biserial r: (T+ - T-) / (T+ + T-), zero differences dropped
    (as in scipy's default Wilcoxon)."""
    d = np.asarray(d, float)
    d = d[d != 0]
    rk = rankdata(np.abs(d))
    return float((rk[d > 0].sum() - rk[d < 0].sum()) / rk.sum())


def _boot_ci(stat, samples, paired, level=0.95):
    from scipy.stats import bootstrap
    rng = np.random.default_rng(BOOT_SEED)
    res = bootstrap(samples, stat, paired=paired, vectorized=False, n_resamples=BOOT_N,
                    confidence_level=level, method='percentile', random_state=rng)
    return float(res.confidence_interval.low), float(res.confidence_interval.high)


def rank_biserial(x, y):
    """Rank-biserial r (x vs y, independent groups) with a percentile bootstrap 95% CI
    (groups resampled separately). Positive = x tends to exceed y."""
    x, y = np.asarray(x, float), np.asarray(y, float)
    return (_rrb_ind(x, y), *_boot_ci(_rrb_ind, (x, y), paired=False))


def rank_biserial_paired(x, y=None):
    """Matched-pairs rank-biserial r for y - x (or for x itself vs 0 when y is None)
    with a percentile bootstrap 95% CI (animals resampled)."""
    d = np.asarray(x, float) if y is None else np.asarray(y, float) - np.asarray(x, float)
    r = _rrb_paired(d)
    if np.all(d > 0) or np.all(d < 0):          # every animal changed in the same direction
        return r, r, r                           # bootstrap is degenerate; CI collapses to +/-1
    return (r, *_boot_ci(_rrb_paired, (d,), paired=False))
