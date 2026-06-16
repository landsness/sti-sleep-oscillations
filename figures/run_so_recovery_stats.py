"""
Compute acute -> week-1 and baseline -> week-1 Wilcoxon signed-rank tests
for each SO power channel, using n=23 (DBSI_M02 and DBSI_M04 excluded).
Analysis only — no figures modified.
"""
import sys
from pathlib import Path
import numpy as np
from scipy.stats import wilcoxon

sys.path.insert(0, str(Path(__file__).parent))
from figure_style import load_data

DATA_FILE = Path(__file__).parent.parent / 'WFCI AI Summary.xlsx'
data = load_data(DATA_FILE)
bl, ac, wk_clean = data['bl'], data['ac'], data['wk_clean']

# Confirm n=23 subset and list animal IDs
idx23 = wk_clean.index.tolist()
print(f"n = {len(idx23)} animals in week-1 subset (DBSI_M02 and DBSI_M04 excluded):")
print("  " + ", ".join(sorted(idx23)))
print()

# Align baseline and acute to the same n=23 animals
bl23 = bl.loc[idx23]
ac23 = ac.loc[idx23]

channels = [
    ('ROI1 Ipsilateral',   'so_power_roi1_ipsi'),
    ('ROI1 Contralateral', 'so_power_roi1_contra'),
    ('ROI2 Ipsilateral',   'so_power_roi2_ipsi'),
    ('ROI2 Contralateral', 'so_power_roi2_contra'),
]

print(f"{'Channel':<22}  {'Acute %BL':>14}  {'Wk1 %BL':>14}  {'Dir':>5}  "
      f"{'Ac->Wk1 W':>10}  {'Ac->Wk1 p':>12}  {'BL->Wk1 W':>10}  {'BL->Wk1 p':>12}")
print("-" * 110)

results = []
for label, col in channels:
    bl_vals  = bl23[col].values
    ac_pct   = ac23[col].values  / bl_vals * 100
    wk_pct   = wk_clean[col].values / bl_vals * 100

    ac_mean, ac_sd = ac_pct.mean(), ac_pct.std(ddof=1)
    wk_mean, wk_sd = wk_pct.mean(), wk_pct.std(ddof=1)

    # Acute -> Week 1 (paired, two-tailed)
    w_ac_wk, p_ac_wk = wilcoxon(ac_pct, wk_pct, alternative='two-sided')

    # Baseline -> Week 1 (i.e., is week-1 still below 100?)
    # Compare wk_pct to a constant 100 — equivalent to wilcoxon(wk_pct - 100)
    w_bl_wk, p_bl_wk = wilcoxon(wk_pct - 100, alternative='two-sided')

    direction = 'UP' if wk_mean > ac_mean else 'DOWN'

    print(f"{label:<22}  {ac_mean:6.1f} ± {ac_sd:5.1f}  {wk_mean:6.1f} ± {wk_sd:5.1f}  "
          f"{direction:>5}  {w_ac_wk:>10.0f}  {p_ac_wk:>12.4e}  {w_bl_wk:>10.0f}  {p_bl_wk:>12.4e}")

    results.append(dict(channel=label, ac_mean=ac_mean, ac_sd=ac_sd,
                        wk_mean=wk_mean, wk_sd=wk_sd, direction=direction,
                        w_ac_wk=w_ac_wk, p_ac_wk=p_ac_wk,
                        w_bl_wk=w_bl_wk, p_bl_wk=p_bl_wk))

print()
print("Notes:")
print("  - All % baseline computed per-animal (channel_value / animal_baseline * 100)")
print("  - Acute -> Wk1: paired Wilcoxon on % baseline values, two-tailed")
print("  - BL -> Wk1: Wilcoxon of (wk_pct - 100), two-tailed; tests whether Wk1 != 100%")
print("  - No multiple-comparison correction (consistent with existing pipeline)")
print(f"  - n = {len(idx23)} for all tests")
