"""
Step 1 — Remaining reviewer-requested analyses for STI_SO revision.
Runs BEFORE FDR correction (these define the correlation family).

Uses the manuscript's own definitions (from figure_style.load_data / make_figure_5):
  LI_roi = (ipsi - contra)/(ipsi + contra)
  Fig 5  = baseline LI_roi1 vs week-1 behavior_score, Spearman, n=25
  partial = rank-residual partial Spearman controlling for stroke_size
  %baseline normalization = (tp / baseline) * 100
  week-1 SO cohort = n=23 (drop DBSI_M02, DBSI_M04); behavior stays n=25
"""
import numpy as np, pandas as pd
from scipy import stats
from scipy.stats import spearmanr, pearsonr, rankdata, linregress, wilcoxon
import statsmodels.formula.api as smf
import statsmodels.api as sm
import json
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
sys.path.insert(0, str(Path(__file__).resolve().parent))
import config  # noqa: E402
from common import corr_ci  # noqa: E402

DATA = config.XLSX
OUTLIERS = ['DBSI_M02', 'DBSI_M04']

df = pd.read_excel(DATA, sheet_name='AI Summary ')
df.columns = df.columns.str.strip()
for c in df.select_dtypes('object').columns:
    df[c] = df[c].str.strip()

bl = df[df.time_point == 'baseline'].set_index('mouse_id').copy()
ac = df[df.time_point == '24 hours'].set_index('mouse_id').copy()
wk = df[df.time_point == '1 week'].set_index('mouse_id').copy()

for tp in (bl, ac, wk):
    for roi in ('roi1', 'roi2'):
        i, c = tp[f'so_power_{roi}_ipsi'], tp[f'so_power_{roi}_contra']
        tp[f'LI_{roi}'] = (i - c) / (i + c)
    for col in [c for c in tp.columns if c.startswith('so_power_')]:
        tp[col + '_pct'] = tp[col] / bl.loc[tp.index, col] * 100

wk_clean = wk.drop([o for o in OUTLIERS if o in wk.index])
shared = bl.index.intersection(wk.index)
bl, wk = bl.loc[shared], wk.loc[shared]
ac = ac.loc[bl.index.intersection(ac.index)]

def partial_spear(x, y, z):
    x, y, z = map(lambda a: np.asarray(a, float), (x, y, z))
    m = ~(np.isnan(x) | np.isnan(y) | np.isnan(z))
    rx, ry, rz = rankdata(x[m]), rankdata(y[m]), rankdata(z[m])
    resid = lambda a, b: a - (linregress(b, a).slope * b + linregress(b, a).intercept)
    r, p = pearsonr(resid(rx, rz), resid(ry, rz))
    return r, p, int(m.sum())

def spear(x, y):
    x, y = np.asarray(x, float), np.asarray(y, float)
    m = ~(np.isnan(x) | np.isnan(y))
    r, p = spearmanr(x[m], y[m])
    return r, p, int(m.sum())

R = {}
line = lambda *a: print(*a)

# ============================================================
# 0. VALIDATION — reproduce known manuscript numbers
# ============================================================
line("="*70); line("VALIDATION — reproduce known manuscript numbers"); line("="*70)
r, p, n = spear(bl['LI_roi1'], wk['behavior_score'])
pr, pp, pn = partial_spear(bl['LI_roi1'], wk['behavior_score'], bl['stroke_size'])
line(f"Fig5A baseline LI_roi1 vs wk1 behavior:  rho={r:.3f} p={p:.3f} n={n}  [target rho=-0.518 p=0.008]")
line(f"Fig5A partial (ctrl infarct):            rho={pr:.3f} p={pp:.3f} n={pn} [target -0.446 p=0.026]")
rb, pb, nb = spear(bl['LI_roi1'], bl['stroke_size'])
line(f"Fig5B baseline LI_roi1 vs infarct vol:   rho={rb:.3f} p={pb:.3f} n={nb}")
# Table 1 check: wk1 SO (%baseline) vs wk1 behavior, lesion ipsi, n=23
wkc = wk_clean
r1, p1, n1 = spear(wkc['so_power_roi1_ipsi_pct'], wkc['behavior_score'])
line(f"Table1 wk1 lesion-ipsi SO(%bl) vs wk1 behavior (n=23): rho={r1:.3f} p={p1:.3f}  [target rho=-0.019]")
R['validation'] = dict(fig5a=(r,p,n), fig5a_partial=(pr,pp,pn), fig5b=(rb,pb,nb), table1_les_ipsi=(r1,p1,n1))

# ============================================================
# R1.13 — infarct volume x acute (24h) SO power
# ============================================================
line("\n" + "="*70); line("R1.13 — infarct volume (stroke_size) x acute 24h SO power"); line("="*70)
acw = ac.loc[bl.index]  # align, n=25
r113 = {}
for label, col in [("lesional ipsi (roi1)", 'so_power_roi1_ipsi'),
                   ("perilesional ipsi (roi2)", 'so_power_roi2_ipsi'),
                   ("lesional contra (roi1)", 'so_power_roi1_contra'),
                   ("perilesional contra (roi2)", 'so_power_roi2_contra')]:
    r_raw, p_raw, n_ = spear(bl['stroke_size'], acw[col])
    r_pct, p_pct, _  = spear(bl['stroke_size'], acw[col + '_pct'])
    pr_r, pr_p = pearsonr(bl['stroke_size'].values, acw[col].values)
    line(f"{label:28s}  raw:  Spearman rho={r_raw:+.3f} p={p_raw:.4f} | Pearson r={pr_r:+.3f} p={pr_p:.4f} (n={n_})")
    line(f"{'':28s}  %bl:  Spearman rho={r_pct:+.3f} p={p_pct:.4f}")
    r113[col] = dict(raw_spear=(r_raw,p_raw), raw_pearson=(pr_r,pr_p), pct_spear=(r_pct,p_pct), n=n_)
R['R1.13'] = r113

# ============================================================
# R1.17 — baseline LI vs baseline forelimb asymmetry (paw preference)
# ============================================================
line("\n" + "="*70); line("R1.17 — baseline LI_roi1 vs baseline forelimb asymmetry"); line("="*70)
r_pp, p_pp, n_pp = spear(bl['LI_roi1'], bl['behavior_score'])
pr_pp, pp_pp = pearsonr(bl['LI_roi1'].values, bl['behavior_score'].values)
line(f"baseline LI_roi1 vs baseline behavior_score: Spearman rho={r_pp:+.3f} p={p_pp:.4f} | Pearson r={pr_pp:+.3f} p={pp_pp:.4f} (n={n_pp})")
# Also roi2
r_pp2, p_pp2, _ = spear(bl['LI_roi2'], bl['behavior_score'])
line(f"baseline LI_roi2 vs baseline behavior_score: Spearman rho={r_pp2:+.3f} p={p_pp2:.4f}")

line("\n-- Partial: baseline LI_roi1 -> wk1 behavior, controlling for baseline forelimb asymmetry --")
pr1, pp1, pn1 = partial_spear(bl['LI_roi1'], wk['behavior_score'], bl['behavior_score'])
ci1_lo, ci1_hi = corr_ci(pr1, pn1, kind='spearman', n_covariates=1)
line(f"partial rho (ctrl baseline asymmetry) = {pr1:+.3f} p={pp1:.4f} n={pn1}  95% CI [{ci1_lo:+.3f}, {ci1_hi:+.3f}]")

line("\n-- wk1 outcome expressed as CHANGE-FROM-BASELINE (wk1 - baseline behavior) --")
chg = (wk['behavior_score'] - bl['behavior_score'])
rc, pc, nc = spear(bl['LI_roi1'], chg)
line(f"baseline LI_roi1 vs change (wk1-bl):        Spearman rho={rc:+.3f} p={pc:.4f} n={nc}")
# partial change controlling infarct
prc, ppc, pnc = partial_spear(bl['LI_roi1'], chg, bl['stroke_size'])
line(f"  partial (ctrl infarct):                   rho={prc:+.3f} p={ppc:.4f}")
R['R1.17'] = dict(li_vs_baseasym=(r_pp,p_pp,n_pp), li_vs_baseasym_pearson=(pr_pp,pp_pp),
                  li2_vs_baseasym=(r_pp2,p_pp2),
                  partial_ctrl_baseasym=(pr1,pp1,pn1),
                  partial_ctrl_baseasym_ci=(ci1_lo,ci1_hi),
                  li_vs_change=(rc,pc,nc), li_vs_change_partial_infarct=(prc,ppc))

# ============================================================
# R2.1 — ANCOVA: does STI add beyond continuous infarct volume?
# ============================================================
line("\n" + "="*70); line("R2.1 — ANCOVA outcome ~ infarct_volume + STI (STI-group outcomes)"); line("="*70)
# collinearity of STI vs infarct volume
mv = bl.groupby('secondary_thalamic_injury')['stroke_size'].agg(['mean','std','count'])
line("Infarct volume by STI status (mm^3):"); line(mv.to_string())
from scipy.stats import mannwhitneyu, pointbiserialr
rpb, ppb = pointbiserialr(bl['secondary_thalamic_injury'], bl['stroke_size'])
line(f"point-biserial STI vs infarct volume: r={rpb:.3f} p={ppb:.2e} (fold-diff means = {mv['mean'][1]/mv['mean'][0]:.1f}x)")

def ancova(name, outcome_series):
    d = pd.DataFrame({'y': outcome_series.values,
                      'infarct': bl['stroke_size'].values,
                      'STI': bl['secondary_thalamic_injury'].values}).dropna()
    line(f"\n--- Outcome: {name} (n={len(d)}) ---")
    # Model 1: infarct only ; Model 2: infarct + STI
    m1 = smf.ols('y ~ infarct', data=d).fit()
    m2 = smf.ols('y ~ infarct + STI', data=d).fit()
    # incremental F for STI beyond infarct
    from statsmodels.stats.anova import anova_lm
    comp = anova_lm(m1, m2)
    dR2 = m2.rsquared - m1.rsquared
    sti_coef = m2.params['STI']; sti_p = m2.pvalues['STI']
    inf_p_full = m2.pvalues['infarct']
    line(f"  infarct-only R2={m1.rsquared:.3f}; +STI R2={m2.rsquared:.3f} (dR2={dR2:.3f})")
    line(f"  STI partial coef={sti_coef:.4g}, p={sti_p:.4f}; infarct p in full model={inf_p_full:.4f}")
    line(f"  incremental F(STI | infarct): F={comp['F'][1]:.3f}, p={comp['Pr(>F)'][1]:.4f}")
    # VIF
    from statsmodels.stats.outliers_influence import variance_inflation_factor
    X = sm.add_constant(d[['infarct','STI']])
    vif = variance_inflation_factor(X.values, 1)
    line(f"  VIF(infarct|STI)={vif:.2f}")
    return dict(n=len(d), r2_infarct=m1.rsquared, r2_full=m2.rsquared, dR2=dR2,
                sti_coef=sti_coef, sti_p=sti_p, infarct_p_full=inf_p_full,
                incF=comp['F'][1], incF_p=comp['Pr(>F)'][1], vif=vif)

acw = ac.loc[bl.index]
r21 = {}
r21['acute_les_ipsi_raw']  = ancova('acute lesional-ipsi SO power (raw)', acw['so_power_roi1_ipsi'])
r21['acute_les_ipsi_pct']  = ancova('acute lesional-ipsi SO suppression (%baseline)', acw['so_power_roi1_ipsi_pct'])
r21['acute_peri_ipsi_pct'] = ancova('acute perilesional-ipsi SO suppression (%baseline)', acw['so_power_roi2_ipsi_pct'])
r21['behavior_24h']        = ancova('24h forelimb asymmetry', acw['behavior_score'])
r21['behavior_wk1']        = ancova('wk1 forelimb asymmetry', wk['behavior_score'])
R['R2.1'] = r21

# ============================================================
# R1.11 — baseline vs 24h behavior (paired Wilcoxon)  [NOT in FDR family]
# ============================================================
line("\n" + "="*70); line("R1.11 — baseline vs 24h forelimb asymmetry (paired Wilcoxon)"); line("="*70)
b = bl['behavior_score']; a = acw['behavior_score']
pair = pd.DataFrame({'bl': b, 'ac': a}).dropna()
W, pW = wilcoxon(pair['bl'], pair['ac'])
line(f"n pairs={len(pair)}; median baseline={pair['bl'].median():.2f}, median 24h={pair['ac'].median():.2f}")
line(f"Wilcoxon signed-rank W={W:.1f}, p={pW:.4f}")
# also |asymmetry| in case magnitude is the deficit measure
Wabs, pWabs = wilcoxon(pair['bl'].abs(), pair['ac'].abs())
line(f"(|asymmetry|) W={Wabs:.1f}, p={pWabs:.4f}; median |bl|={pair['bl'].abs().median():.2f}, |24h|={pair['ac'].abs().median():.2f}")
R['R1.11'] = dict(n=len(pair), W=float(W), p=float(pW), med_bl=float(pair['bl'].median()),
                  med_ac=float(pair['ac'].median()), Wabs=float(Wabs), pabs=float(pWabs))

config.STATS.mkdir(parents=True, exist_ok=True)
with open(config.STATS / 'step1_results.json', 'w') as f:
    json.dump(R, f, indent=2, default=str)
line("\n[saved results json]")
