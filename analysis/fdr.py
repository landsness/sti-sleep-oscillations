"""fdr.py -- multiple-comparison control for the correlation analyses (Reviewer 2, point 2).

Family definition (fixed before any correction was computed; see Methods):
  * CONFIRMATORY family (24 Spearman / partial-Spearman correlations) addressing
    the study's two primary questions -- (1) acute SO suppression vs injury and
    concurrent deficit, (2) SO vs week-1 outcome. p-values are adjusted with the
    Benjamini-Hochberg procedure and reported as q.
  Each rho is reported with its 95% CI (Fisher z, Bonett-Wright SE; Reviewer 2, minor 2).
  * EXPLORATORY family (10 correlations) -- pre-stroke (baseline) SO power and
    laterality index vs outcome, labelled exploratory at Reviewer 2's request;
    reported with uncorrected p. BH q within this block, and across all 34
    correlations, are provided for transparency in Supplementary Table 2.
Group / paired comparisons (Mann-Whitney U, Wilcoxon) and bilateral-covariation
tests are planned contrasts and are not part of either family.

Figure scripts call `q_value(key, p_check)` so every q printed on a figure comes
from this single table, and the figure's own p is asserted to match.

Output: results/tables/SupplTable2_FDR.csv
"""
import sys
from functools import lru_cache
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import spearmanr
from statsmodels.stats.multitest import multipletests

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
sys.path.insert(0, str(Path(__file__).resolve().parent))
import config  # noqa: E402
from common import load_data, partial_spearman, corr_ci  # noqa: E402

CH = {'lesion ipsi': 'so_power_roi1_ipsi', 'peri ipsi': 'so_power_roi2_ipsi',
      'lesion contra': 'so_power_roi1_contra', 'peri contra': 'so_power_roi2_contra'}
ROI = {'lesion': 'roi1', 'peri': 'roi2'}


def _tests():
    """Yield (family, key, description, figure/text location, n, rho, p)."""
    d = load_data()
    bl = d['bl']
    ac = d['ac'].loc[bl.index]
    wk = d['wk'].loc[bl.index]
    wkc = d['wk_clean']
    inf = bl['stroke_size']
    sti = bl['secondary_thalamic_injury']
    S = lambda x, y: spearmanr(np.asarray(x, float), np.asarray(y, float))

    C = 'confirmatory'
    # -- infarct vs behaviour (Suppl Fig 2)
    yield (C, 'infarct_vs_acute_behavior', 'Infarct volume vs acute behavior', 'Suppl Fig 2A; §3.1', 25, *S(inf, ac.behavior_score))
    yield (C, 'infarct_vs_wk1_behavior', 'Infarct volume vs week-1 behavior', 'Suppl Fig 2B; §3.1', 25, *S(inf, wk.behavior_score))
    # -- infarct vs acute SO power, raw (Suppl Fig 3)
    for lab, col in CH.items():
        yield (C, f'infarct_vs_acute_so_{lab}', f'Infarct volume vs acute SO power, {lab}', 'Suppl Fig 3; §3.2', 25, *S(inf, ac[col]))
    # -- acute SO power vs acute behaviour (Fig 2F-I)
    for lab, col in CH.items():
        yield (C, f'acute_so_vs_acute_behavior_{lab}', f'Acute SO power vs acute behavior, {lab}', 'Fig 2F-I; §3.2; Abstract', 25, *S(ac[col], ac.behavior_score))
    # -- partial (controlling infarct volume)
    for lab in ('lesion ipsi', 'peri ipsi'):
        yield (C, f'partial_acute_so_vs_acute_behavior_{lab}', f'Acute SO power vs acute behavior | infarct, {lab}', '§3.2; Abstract', 25,
               *partial_spearman(ac[CH[lab]], ac.behavior_score, inf))
    # -- STI subgroups (Suppl Fig 4)
    for g, m in (('STI+', sti == 1), ('STI-', sti == 0)):
        yield (C, f'sti_subgroup_{g}_acute_so_vs_acute_behavior', f'{g} subgroup: acute lesion-ipsi SO vs acute behavior', 'Suppl Fig 4; §3.2', int(m.sum()),
               *S(ac.loc[m, CH['lesion ipsi']], ac.loc[m, 'behavior_score']))
    # -- acute LI (Suppl Fig 5)
    for lab, r in ROI.items():
        yield (C, f'acute_li_vs_acute_behavior_{lab}', f'Acute LI vs acute behavior, {lab}', 'Suppl Fig 5A-B; §3.2', 25, *S(ac[f'LI_{r}'], ac.behavior_score))
    for lab, r in ROI.items():
        yield (C, f'acute_li_vs_wk1_behavior_{lab}', f'Acute LI vs week-1 behavior, {lab}', 'Suppl Fig 5C-D; §3.2', 25, *S(ac[f'LI_{r}'], wk.behavior_score))
    for lab, r in ROI.items():
        yield (C, f'acute_li_vs_infarct_{lab}', f'Acute LI vs infarct volume, {lab}', 'Suppl Fig 5E-F; §3.2', 25, *S(ac[f'LI_{r}'], inf))
    # -- week-1 SO (% baseline) vs week-1 behaviour, n = 23 (Suppl Fig 6; formerly Table 1)
    for lab, col in CH.items():
        pct = wkc[col] / bl.loc[wkc.index, col] * 100
        yield (C, f'wk1_so_vs_wk1_behavior_{lab}', f'Week-1 SO (% baseline) vs week-1 behavior, {lab}', 'Suppl Fig 6; §3.3; Abstract', 23, *S(pct, wkc.behavior_score))

    E = 'exploratory'
    for lab, col in CH.items():
        yield (E, f'baseline_so_vs_wk1_behavior_{lab}', f'Baseline SO power vs week-1 behavior, {lab}', 'Suppl Fig 7; §3.4', 25, *S(bl[col], wk.behavior_score))
    yield (E, 'baseline_li_vs_wk1_behavior_lesion', 'Baseline LI vs week-1 behavior, lesion', 'Fig 5A; §3.4; Abstract', 25, *S(bl.LI_roi1, wk.behavior_score))
    yield (E, 'partial_baseline_li_vs_wk1_behavior_lesion', 'Baseline LI vs week-1 behavior | infarct, lesion', 'Fig 5A; §3.4; Abstract', 25,
           *partial_spearman(bl.LI_roi1, wk.behavior_score, inf))
    yield (E, 'baseline_li_vs_infarct_lesion', 'Baseline LI vs infarct volume, lesion', 'Fig 5B; §3.4', 25, *S(bl.LI_roi1, inf))
    yield (E, 'baseline_li_vs_wk1_behavior_peri', 'Baseline LI vs week-1 behavior, perilesional', 'Suppl Fig 8; §3.4', 25, *S(bl.LI_roi2, wk.behavior_score))
    yield (E, 'partial_baseline_li_vs_wk1_behavior_peri', 'Baseline LI vs week-1 behavior | infarct, perilesional', '§3.4', 25,
           *partial_spearman(bl.LI_roi2, wk.behavior_score, inf))
    yield (E, 'baseline_li_vs_baseline_asymmetry_lesion', 'Baseline LI vs baseline forelimb asymmetry (paw preference, R1.17)', '§3.4', 25,
           *S(bl.LI_roi1, bl.behavior_score))


@lru_cache(maxsize=1)
def table():
    df = pd.DataFrame(list(_tests()), columns=['family', 'key', 'test', 'reported_in', 'n', 'rho', 'p'])
    assert df.key.is_unique
    assert (df.family == 'confirmatory').sum() == 24 and (df.family == 'exploratory').sum() == 10
    # 95% CI for rho: Fisher z with Bonett-Wright SE; one covariate for the partial correlations
    ci = [corr_ci(r.rho, r.n, n_covariates=int(r.key.startswith('partial'))) for r in df.itertuples()]
    df.insert(df.columns.get_loc('rho') + 1, 'ci95_lo', [c[0] for c in ci])
    df.insert(df.columns.get_loc('rho') + 2, 'ci95_hi', [c[1] for c in ci])
    bh = lambda p: multipletests(p, method='fdr_bh')[1]
    df['q'] = np.nan                                    # BH within the confirmatory family (reported)
    m = df.family == 'confirmatory'
    df.loc[m, 'q'] = bh(df.loc[m, 'p'])
    df['q_within_exploratory'] = np.nan
    df.loc[~m, 'q_within_exploratory'] = bh(df.loc[~m, 'p'])
    df['q_all_34'] = bh(df['p'])
    return df


def q_value(key, p_check=None, tol=1e-9):
    """BH q for a confirmatory test (NaN for exploratory). If p_check is given,
    assert it equals the family p for that key (guards figure/table drift)."""
    row = table().set_index('key').loc[key]
    if p_check is not None and abs(row.p - p_check) > tol:
        raise AssertionError(f'{key}: figure p={p_check} != family p={row.p}')
    return row.q


def main():
    df = table()
    config.TABLES.mkdir(parents=True, exist_ok=True)
    out = config.TABLES / 'SupplTable2_FDR.csv'
    df.to_csv(out, index=False)
    pd.set_option('display.width', 200, 'display.max_colwidth', 70)
    print(df[['family', 'test', 'n', 'rho', 'ci95_lo', 'ci95_hi', 'p', 'q', 'q_within_exploratory', 'q_all_34']].round(4).to_string(index=False))
    print(f'\nSaved {out}')


if __name__ == '__main__':
    main()
