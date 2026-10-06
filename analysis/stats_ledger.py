"""stats_ledger.py -- Recompute every inferential test reported in the manuscript
(reported values as of _092926v4.docx) from WFCI_AI_Summary.xlsx using the pipeline's own
definitions (common.load_data), compare to the reported values, and attach the effect size
with its 95% CI (Reviewer 2, minor point 2):
  * Spearman rho / partial rho: Fisher z, Bonett-Wright SE (common.corr_ci)
  * Pearson r (bilateral covariation): Fisher z
  * Mann-Whitney U: rank-biserial r, percentile bootstrap CI (common.rank_biserial)
  * Wilcoxon signed-rank: matched-pairs rank-biserial r, bootstrap CI (common.rank_biserial_paired)
Effect sizes are captured automatically from the test call on the same line (see the
wrappers below), so each row's effect size is computed on exactly the data of its test."""
import sys, numpy as np, pandas as pd
from pathlib import Path
HERE = Path(__file__).resolve().parent
sys.path[:0] = [str(HERE), str(HERE.parent)]
import config  # noqa: E402
from common import load_data, corr_ci, rank_biserial, rank_biserial_paired  # noqa: E402
from scipy.stats import rankdata, linregress, norm
from scipy import stats as _st

# -- thin wrappers: run the test and remember its effect size + 95% CI for the next add()
_ES = {}
def spearmanr(x, y):
    r, p = _st.spearmanr(x, y); _ES['es'] = ('rho', r, *corr_ci(r, len(x))); return r, p
def pearsonr(x, y):
    r, p = _st.pearsonr(x, y); _ES['es'] = ('r', r, *corr_ci(r, len(x), kind='pearson')); return r, p
def mannwhitneyu(x, y):
    res = _st.mannwhitneyu(x, y); _ES['es'] = ('r_rb', *rank_biserial(x, y)); return res
def wilcoxon(x, y=None):
    res = _st.wilcoxon(x) if y is None else _st.wilcoxon(x, y)
    _ES['es'] = ('r_rb_paired', *rank_biserial_paired(x, y)); return res
import statsmodels.formula.api as smf
from statsmodels.stats.anova import anova_lm

D = load_data(config.XLSX)
bl, ac, wk, wkc = D['bl'], D['ac'].loc[D['bl'].index], D['wk'].loc[D['bl'].index], D['wk_clean']
idx23 = wkc.index
sti = bl['secondary_thalamic_injury']
P, N = sti == 1, sti == 0
inf = bl['stroke_size']
CH = {'L-ipsi': 'so_power_roi1_ipsi', 'P-ipsi': 'so_power_roi2_ipsi',
      'L-contra': 'so_power_roi1_contra', 'P-contra': 'so_power_roi2_contra'}

def pspear(x, y, z):
    x, y, z = [np.asarray(v, float) for v in (x, y, z)]
    rx, ry, rz = rankdata(x), rankdata(y), rankdata(z)
    res = lambda a, b: a - (linregress(b, a).slope * b + linregress(b, a).intercept)
    r, p = _st.pearsonr(res(rx, rz), res(ry, rz))
    _ES['es'] = ('partial_rho', r, *corr_ci(r, len(x), n_covariates=1)); return r, p

def fisher_dep_naive(r1, r2, n1, n2):  # independent-samples Fisher z (as the manuscript describes)
    z = (np.arctanh(r1) - np.arctanh(r2)) / np.sqrt(1/(n1-3) + 1/(n2-3))
    return z, 2*norm.sf(abs(z))

def steiger_repo(r12, r34, r13, r14, r23, r24, n):
    """Exact copy of make_suppl_table_sensitivity.steiger1980 (the implementation behind the manuscript numbers)."""
    z12, z34 = np.arctanh(r12), np.arctanh(r34)
    rbar = (r12 + r34) / 2
    num = (0.5 * r12 * r34 * (r13**2 + r14**2 + r23**2 + r24**2)
           + r13 * r24 + r14 * r23
           - (r12 * r13 * r14 + r12 * r23 * r24 + r13 * r23 * r34 + r14 * r24 * r34))
    c = num / ((1 - rbar**2)**2)
    Z = (z12 - z34) * np.sqrt((n - 3) / (2 - 2 * c))
    return Z, 2 * norm.sf(abs(Z))

def steiger1980(r12, r34, r13, r14, r23, r24, n):
    """Steiger (1980) eq. with pooled rbar substituted throughout psi (textbook form)."""
    z12, z34 = np.arctanh(r12), np.arctanh(r34)
    rbar = (r12 + r34) / 2
    psi = 0.5*((r13 - rbar*r23)*(r24 - r23*r34) + (r14 - r13*r34)*(r23 - rbar*r13)
               + (r13 - r14*rbar)*(r24 - rbar*r14) + (r14 - rbar*r24)*(r23 - r24*rbar))
    c = psi / (1 - rbar**2)**2
    Z = (z12 - z34) * np.sqrt(n - 3) / np.sqrt(2 - 2*c)
    return Z, 2*norm.sf(abs(Z))

rows = []
def add(id_, sec, test, family, stat, p, rep_stat, rep_p, script, n, note=''):
    es = _ES.pop('es', None)
    if family in ('Covariation change', 'Model comparison'):
        es = None                                  # effect size is the correlations / delta-R2 themselves
    es_type, es_val, lo, hi = es if es else ('', np.nan, np.nan, np.nan)
    rows.append(dict(id=id_, section=sec, test=test, family=family, stat=round(float(stat), 3),
                     p_computed=float(p), reported_stat=rep_stat, reported_p=rep_p, n=n,
                     effect_size_type=es_type, effect_size=es_val, ci95_lo=lo, ci95_hi=hi,
                     source_script=script, note=note))

pct = lambda df, c, ix: (df.loc[ix, c] / bl.loc[ix, c] * 100)

# ---------------- §3.1 ----------------
u, p = mannwhitneyu(inf[P], inf[N]); add('T01', '3.1', 'Infarct volume STI+ vs STI-', 'Group (STI)', u, p, 'U=151', '<0.001', 'make_figure_1.py', 25)
r, p = spearmanr(inf, ac['behavior_score']); add('T02', '3.1', 'Infarct vs acute behavior', 'Corr: structure-behavior', r, p, 'ρ=0.636', '<0.001', 'make_suppl_s1.py', 25)
r, p = spearmanr(inf, wk['behavior_score']); add('T03', '3.1', 'Infarct vs wk1 behavior', 'Corr: structure-behavior', r, p, 'ρ=0.502', '0.011', 'make_suppl_s1.py', 25)
u, p = mannwhitneyu(bl.loc[P, 'behavior_score'], bl.loc[N, 'behavior_score']); add('T04', '3.1', 'Baseline behavior STI+ vs STI-', 'Group (STI)', u, p, 'U=82', '0.805', 'make_figure_1.py', 25)
w, p = wilcoxon(bl['behavior_score'], ac['behavior_score']); add('T05', '3.1', 'Behavior baseline vs 24h (R1.11)', 'Manipulation check', w, p, 'W=21', '<0.001', 'step1_remaining_analyses.py', 25)
u, p = mannwhitneyu(ac.loc[P, 'behavior_score'], ac.loc[N, 'behavior_score']); add('T06', '3.1', 'Acute behavior STI+ vs STI-', 'Group (STI)', u, p, 'U=106', '0.119', 'make_figure_1.py', 25)

# ---------------- §3.2 acute suppression ----------------
for k, c in CH.items():
    w, p = wilcoxon(pct(ac, c, bl.index) - 100)
    add(f'T07{k}', '3.2', f'Acute vs baseline, {k} (%BL vs 100)', 'Within-animal SO', w, p, '—', '<0.001 (all)', 'NONE', 25, 'text says "all Wilcoxon p<0.001"; no script computes it')
w, p = wilcoxon(pct(ac, CH['L-ipsi'], bl.index), pct(ac, CH['L-contra'], bl.index)); add('T08', '3.2', 'Acute ipsi vs contra, lesion', 'Within-animal SO', w, p, 'W=3', '<0.001', 'make_figure_2.py (panel A)', 25)
w, p = wilcoxon(pct(ac, CH['P-ipsi'], bl.index), pct(ac, CH['P-contra'], bl.index)); add('T09', '3.2', 'Acute ipsi vs contra, peri', 'Within-animal SO', w, p, 'W=8', '<0.001', 'make_figure_2.py (panel A)', 25)
w, p = wilcoxon(pct(ac, CH['L-ipsi'], bl.index), pct(ac, CH['P-ipsi'], bl.index)); add('T10', '3.2', 'Acute lesion vs peri ipsi', 'Within-animal SO', w, p, 'W=1', '<0.001', 'make_figure_2.py (panel A)', 25,
    f"means {pct(ac, CH['L-ipsi'], bl.index).mean():.1f}% vs {pct(ac, CH['P-ipsi'], bl.index).mean():.1f}% (text 17.2 vs 30.8)")
for k, rep in zip(CH, [('U=38', '0.035'), ('U=45', '0.085'), ('U=78', '0.978'), ('U=93', '0.396')]):
    u, p = mannwhitneyu(pct(ac, CH[k], bl.index)[P], pct(ac, CH[k], bl.index)[N])
    add(f'T11{k}', '3.2', f'Acute SO %BL STI+ vs STI-, {k}', 'Group (STI)', u, p, rep[0], rep[1], 'make_figure_2.py (panel B)', 25)
for k, rep in zip(CH, [('ρ=−0.809', '<0.001'), ('ρ=−0.742', '<0.001'), ('ρ=0.092', '0.663'), ('ρ=−0.078', '0.712')]):
    r, p = spearmanr(inf, ac[CH[k]]); add(f'T12{k}', '3.2', f'Infarct vs acute raw SO, {k} (R1.13)', 'Corr: structure-SO', r, p, rep[0], rep[1], 'make_suppl_infarct_vs_acuteSO.py', 25)
d = pd.DataFrame({'y': pct(ac, CH['L-ipsi'], bl.index).values, 'infv': inf.values, 'sti': sti.values})
m1, m2 = smf.ols('y~infv', d).fit(), smf.ols('y~infv+sti', d).fit(); cmp_ = anova_lm(m1, m2)
add('T13', '3.2', 'Nested regression: STI beyond infarct (lesion-ipsi %BL) (R2.1)', 'Model comparison', m2.rsquared - m1.rsquared, cmp_['Pr(>F)'][1], 'ΔR²=0.01', '0.566', 'step1_remaining_analyses.py', 25,
    'regression uses %BL; adjacent infarct correlations (T12) use raw power')
for k, rep in zip(CH, [('ρ=−0.718', '<0.001'), ('ρ=−0.600', '0.002'), ('ρ=−0.047', '0.824'), ('ρ=−0.079', '0.707')]):
    r, p = spearmanr(ac[CH[k]], ac['behavior_score']); add(f'T14{k}', '3.2', f'Acute raw SO vs acute behavior, {k}', 'Corr: SO-behavior', r, p, rep[0], rep[1], 'make_figure_2.py', 25)
r, p = pspear(ac[CH['L-ipsi']], ac['behavior_score'], inf); add('T15', '3.2', 'Partial: acute lesion-ipsi SO vs behavior | infarct', 'Corr: SO-behavior (partial)', r, p, 'ρ=−0.448', '0.025', 'NONE', 25, 'headline abstract stat; not computed in any repo script')
r, p = pspear(ac[CH['P-ipsi']], ac['behavior_score'], inf); add('T16', '3.2', 'Partial: acute peri-ipsi SO vs behavior | infarct', 'Corr: SO-behavior (partial)', r, p, 'ρ=−0.248', '0.232', 'NONE', 25)
r, p = spearmanr(ac.loc[P, CH['L-ipsi']], ac.loc[P, 'behavior_score']); add('T17a', '3.2', 'STI+ subgroup: acute lesion-ipsi SO vs behavior', 'Corr: SO-behavior (subgroup)', r, p, 'ρ=−0.666', '0.009', 'make_suppl_s2.py', int(P.sum()))
r, p = spearmanr(ac.loc[N, CH['L-ipsi']], ac.loc[N, 'behavior_score']); add('T17b', '3.2', 'STI- subgroup: acute lesion-ipsi SO vs behavior', 'Corr: SO-behavior (subgroup)', r, p, 'ρ=−0.673', '0.023', 'make_suppl_s2.py', int(N.sum()))
for roi, rep in [('roi1', ('ρ=−0.649', '<0.001')), ('roi2', ('ρ=−0.612', '0.001'))]:
    r, p = spearmanr(ac[f'LI_{roi}'], ac['behavior_score']); add(f'T18{roi}', '3.2', f'Acute LI_{roi} vs acute behavior', 'Corr: LI-behavior', r, p, rep[0], rep[1], 'make_suppl_s3.py', 25)
for roi, rep in [('roi1', ('ρ=−0.482', '0.015')), ('roi2', ('ρ=−0.516', '0.008'))]:
    r, p = spearmanr(ac[f'LI_{roi}'], wk['behavior_score']); add(f'T19{roi}', '3.2', f'Acute LI_{roi} vs wk1 behavior', 'Corr: LI-behavior', r, p, rep[0], rep[1], 'make_suppl_s3.py', 25)
for roi, rep in [('roi1', ('ρ=−0.744', '<0.001')), ('roi2', ('ρ=−0.675', '<0.001'))]:
    r, p = spearmanr(ac[f'LI_{roi}'], inf); add(f'T20{roi}', '3.2', f'Acute LI_{roi} vs infarct', 'Corr: structure-SO', r, p, rep[0], rep[1], 'make_suppl_s3.py', 25)
# bilateral covariation
cov = {}
for tp, df in [('bl', bl), ('ac', ac)]:
    for roi in ['roi1', 'roi2']:
        r, p = pearsonr(df[f'so_power_{roi}_ipsi'], df[f'so_power_{roi}_contra']); cov[(tp, roi)] = r
        rep = {('bl', 'roi1'): ('r=0.865', '<0.001'), ('bl', 'roi2'): ('r=0.922', '<0.001'), ('ac', 'roi1'): ('r=0.213', '0.307'), ('ac', 'roi2'): ('r=0.475', '0.016')}[(tp, roi)]
        add(f'T21{tp}{roi}', '3.2', f'Bilateral covariation {tp} {roi}', 'Covariation', r, p, rep[0], rep[1], 'make_figure_3.py', 25)
for roi, rep in [('roi1', ('z=3.64', '<0.001')), ('roi2', ('z=3.60', '<0.001'))]:
    z, p = fisher_dep_naive(cov[('bl', roi)], cov[('ac', roi)], 25, 25)
    add(f'T22{roi}', '3.2', f'Covariation baseline vs acute, {roi} (Fisher z, independent)', 'Covariation change', z, p, rep[0], rep[1], 'NONE', 25,
        'independent-samples Fisher z applied to the SAME animals; acute->wk1 uses Steiger (dependent). Inconsistent')
    ix = bl.index
    i, c = f'so_power_{roi}_ipsi', f'so_power_{roi}_contra'
    R = lambda a, b: pearsonr(a, b)[0]
    zs, ps = steiger_repo(cov[('bl', roi)], cov[('ac', roi)], R(bl[i], ac[i]), R(bl[i], ac[c]), R(bl[c], ac[i]), R(bl[c], ac[c]), 25)
    add(f'T22{roi}s', '3.2', f'  (same, Steiger dependent test as used for acute->wk1 — for comparison)', 'Covariation change', zs, ps, '—', '—', 'this audit', 25)
for g, m in [('STI+', P), ('STI-', N)]:
    for roi in ['roi1', 'roi2']:
        r, p = pearsonr(ac.loc[m, f'so_power_{roi}_ipsi'], ac.loc[m, f'so_power_{roi}_contra'])
        add(f'T23{g}{roi}', '3.2', f'Acute covariation within {g}, {roi}', 'Covariation (subgroup)', r, p, '—', '>0.16 (all)', 'NONE', int(m.sum()))

# ---------------- §3.3 recovery (n=23) ----------------
for k, rep in zip(CH, [('W=7', '<0.001'), ('W=10', '<0.001'), ('W=81', '0.086'), ('W=106', '0.345')]):
    w, p = wilcoxon(pct(ac, CH[k], idx23), pct(wkc, CH[k], idx23)); add(f'T24{k}', '3.3', f'Acute->wk1 SO, {k}', 'Within-animal SO', w, p, rep[0], rep[1], 'run_so_recovery_stats.py / make_suppl_table_sensitivity.py', 23)
for k in CH:
    w, p = wilcoxon(pct(wkc, CH[k], idx23) - 100); add(f'T25{k}', '3.3', f'Wk1 vs baseline, {k}', 'Within-animal SO', w, p, 'W≤3', '<0.001 (all)', 'run_so_recovery_stats.py', 23)
w, p = wilcoxon(pct(wkc, CH['L-ipsi'], idx23), pct(wkc, CH['P-ipsi'], idx23)); add('T26', '3.3', 'Wk1 lesion vs peri ipsi', 'Within-animal SO', w, p, 'W=24', '<0.001', 'make_suppl_table_sensitivity.py', 23)
w, p = wilcoxon(pct(wkc, CH['L-ipsi'], idx23), pct(wkc, CH['L-contra'], idx23)); add('T27', '3.3', 'Wk1 ipsi vs contra, lesion', 'Within-animal SO', w, p, 'W=52', '0.007', 'make_suppl_table_sensitivity.py', 23)
w, p = wilcoxon(pct(wkc, CH['P-ipsi'], idx23), pct(wkc, CH['P-contra'], idx23)); add('T28', '3.3', 'Wk1 ipsi vs contra, peri', 'Within-animal SO', w, p, 'W=120', '0.601', 'NONE', 23)
for roi, rep in [('roi1', 'r=0.733'), ('roi2', 'r=0.838')]:
    i, c = f'so_power_{roi}_ipsi', f'so_power_{roi}_contra'
    r, p = pearsonr(wkc[i], wkc[c]); add(f'T29{roi}', '3.3', f'Wk1 covariation {roi}', 'Covariation', r, p, rep, '<0.001 (Fig 4B)', 'make_figure_4.py', 23)
    a, w_ = ac.loc[idx23], wkc
    R = lambda x, y: pearsonr(x, y)[0]
    args = (R(a[i], a[c]), R(w_[i], w_[c]), R(a[i], w_[i]), R(a[i], w_[c]), R(a[c], w_[i]), R(a[c], w_[c]), 23)
    zp, pp = steiger1980(*args)
    add(f'T30{roi}s', '3.3', f'  (same, pooled-rbar Steiger form — for comparison)', 'Covariation change', zp, pp, '—', '—', 'this audit', 23)
    z, p = steiger_repo(*args)
    add(f'T30{roi}', '3.3', f'Covariation acute->wk1 Steiger, {roi}', 'Covariation change', z, p, {'roi1': 'z=−2.49', 'roi2': 'z=−2.37'}[roi], {'roi1': '0.013', 'roi2': '0.018'}[roi], 'make_suppl_table_sensitivity.py', 23)
for k, rep in zip(CH, [('ρ=−0.019', '0.932'), ('ρ=−0.050', '0.819'), ('ρ=−0.074', '0.737'), ('ρ=+0.016', '0.943')]):
    r, p = spearmanr(pct(wkc, CH[k], idx23), wkc['behavior_score']); add(f'T31{k}', '3.3', f'Wk1 SO %BL vs wk1 behavior, {k} (Suppl Fig 6)', 'Corr: SO-behavior', r, p, rep[0], rep[1], 'make_suppl_table_sensitivity.py', 23)
Pc, Nc = D['sti_pos_clean'], D['sti_neg_clean']
for k in ['L-ipsi', 'P-ipsi']:
    v = pct(wkc, CH[k], idx23)
    u, p = mannwhitneyu(v[Pc.values], v[Nc.values])
    add(f'T32{k}', '3.3', f'Wk1 SO recovery STI+ vs STI-, {k}', 'Group (STI)', u, p, 'U=69' if k == 'L-ipsi' else '—', '0.878' if k == 'L-ipsi' else '"any region"',
        'make_figure_4.py' if k == 'L-ipsi' else 'NONE', 23, f'means {v[Pc.values].mean():.1f}% vs {v[Nc.values].mean():.1f}% (text: 41.8 vs 38.8)' if k == 'L-ipsi' else 'text claims "any region"; peri not computed')
u, p = mannwhitneyu(wk.loc[P, 'behavior_score'], wk.loc[N, 'behavior_score']); add('T33', '3.3', 'Wk1 behavior STI+ vs STI-', 'Group (STI)', u, p, 'U=124', '0.011', 'make_figure_4.py', 25)

# ---------------- §3.4 baseline ----------------
for k, rep in zip(CH, [('ρ=−0.135', '0.519'), ('ρ=−0.188', '0.369'), ('ρ=0.165', '0.432'), ('ρ=−0.143', '0.495')]):
    r, p = spearmanr(bl[CH[k]], wk['behavior_score']); add(f'T34{k}', '3.4', f'Baseline SO vs wk1 behavior, {k}', 'Corr: SO-behavior', r, p, rep[0], rep[1], 'make_suppl_s4.py', 25)
r, p = spearmanr(bl['LI_roi1'], wk['behavior_score']); add('T35', '3.4', 'Baseline LI lesion vs wk1 behavior', 'Corr: LI-behavior', r, p, 'ρ=−0.518', '0.008', 'make_figure_5.py', 25)
r, p = pspear(bl['LI_roi1'], wk['behavior_score'], inf); add('T36', '3.4', 'Partial: baseline LI lesion vs wk1 | infarct', 'Corr: LI-behavior (partial)', r, p, 'ρ=−0.446', '0.026', 'make_figure_5.py', 25)
r, p = spearmanr(bl['LI_roi1'], inf); add('T37', '3.4', 'Baseline LI lesion vs infarct', 'Corr: structure-SO', r, p, 'ρ=−0.298', '0.147', 'make_figure_5.py', 25)
r, p = spearmanr(bl['LI_roi2'], wk['behavior_score']); add('T38', '3.4', 'Baseline LI peri vs wk1 behavior', 'Corr: LI-behavior', r, p, 'ρ=−0.021', '0.922', 'make_suppl_s5.py', 25)
r, p = pspear(bl['LI_roi2'], wk['behavior_score'], inf); add('T39', '3.4', 'Partial: baseline LI peri vs wk1 | infarct', 'Corr: LI-behavior (partial)', r, p, 'ρ=−0.092', '0.664', 'NONE', 25)
r, p = spearmanr(bl['LI_roi1'], bl['behavior_score']); add('T40', '3.4', 'Baseline LI lesion vs baseline asymmetry (R1.17)', 'Corr: LI-behavior', r, p, 'ρ=−0.011', '0.959', 'step1_remaining_analyses.py', 25)
r, p = pspear(bl['LI_roi1'], wk['behavior_score'], bl['behavior_score']); add('T41', '3.4', 'Partial: baseline LI lesion vs wk1 | baseline asymmetry (R1.17)', 'Corr: LI-behavior (partial)', r, p, 'ρ=−0.542', '0.005', 'step1_remaining_analyses.py', 25)

df = pd.DataFrame(rows)

def parse_p(s):
    s = str(s).replace('(all)', '').replace('(Fig 4B)', '').strip()
    try:
        if s.startswith('<'): return ('lt', float(s[1:]))
        if s.startswith('>'): return ('gt', float(s[1:]))
        return ('eq', float(s))
    except: return (None, None)
def check(r):
    kind, v = parse_p(r.reported_p)
    pc = r.p_computed
    if kind is None: return 'n/a'
    if kind == 'lt': return 'OK' if pc < v else 'MISMATCH'
    if kind == 'gt': return 'OK' if pc > v else 'MISMATCH'
    return 'OK' if abs(pc - v) < 0.0015 or (v >= 0.1 and abs(pc - v) < 0.006) else 'MISMATCH'
df['p_match'] = df.apply(check, axis=1)
FAM = {
 'Corr: SO-behavior':'F1 correlations', 'Corr: SO-behavior (partial)':'F1 correlations', 'Corr: SO-behavior (subgroup)':'F1 correlations',
 'Corr: LI-behavior':'F1 correlations', 'Corr: LI-behavior (partial)':'F1 correlations', 'Corr: structure-SO':'F1 correlations',
 'Corr: structure-behavior':'F1 correlations',
 'Covariation':'F2 covariation', 'Covariation change':'F2 covariation', 'Covariation (subgroup)':'F2 covariation',
 'Group (STI)':'F3 group/paired', 'Within-animal SO':'F3 group/paired',
 'Manipulation check':'none (manipulation check)', 'Model comparison':'none (confound model, R2.1)'}
df['fdr_family_proposed'] = df['family'].map(FAM)
df.loc[df.source_script=='this audit','fdr_family_proposed'] = 'n/a (audit only)'
df.loc[df.id=='T40','fdr_family_proposed'] = 'F1 correlations'
df['trace_status'] = np.where(df.source_script.eq('NONE'), 'UNTRACED', 'traced')
pd.set_option('display.width', 250, 'display.max_colwidth', 60, 'display.max_rows', 200)
config.STATS.mkdir(parents=True, exist_ok=True)
OUT = config.STATS / 'stats_ledger.csv'
df.to_csv(OUT, index=False)
print(df[['id', 'test', 'stat', 'p_computed', 'reported_stat', 'reported_p', 'p_match', 'effect_size_type', 'effect_size', 'ci95_lo', 'ci95_hi']].round(3).to_string())
print('\nNOTES'); print(df.loc[df.note != '', ['id', 'note']].to_string())
