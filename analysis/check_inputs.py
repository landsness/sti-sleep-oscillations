"""check_inputs.py -- verify the hand-curated analysis table against its machine-generated sources.

data/WFCI_AI_Summary.xlsx was assembled by hand. This check confirms that every
value the statistics depend on matches the upstream outputs it was copied from:
  * SO power (4 channels x 3 timepoints x 25 mice = 300 values)
        vs data/ROI_Summary_V2.xlsx  (upstream/matlab/SO_ROI_Quantification_mirrored_v2.m)
  * infarct volume (stroke_size, mm^3; 25 mice)
        vs data/Stroke_Quantification.xlsx  (per-slice ImageJ infarct areas x 40 um)
Exits non-zero if any value differs by more than the tolerance.

Output: results/stats/input_checks.csv
"""
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import config  # noqa: E402
from common import load_data, roi_mouse_id  # noqa: E402

REL_TOL = 1e-4
STROKE_XLSX = config.DATA / 'Stroke_Quantification.xlsx'
TP = {'bl': 'bsl', 'ac': 'acute', 'wk': 'wk1'}


def check_so_power(d):
    roi = pd.read_excel(config.ROI_SUMMARY)
    roi.index = roi['Mouse'].map(roi_mouse_id)
    rows = []
    for key, tp in TP.items():
        df = d[key]
        for mid in df.index:
            for r in ('roi1', 'roi2'):
                for h, H in (('ipsi', 'Ipsi'), ('contra', 'Contra')):
                    a = df.loc[mid, f'so_power_{r}_{h}']
                    b = roi.loc[mid, f'{r.upper()}_{tp}_{H}']
                    rows.append(dict(check='SO power', mouse_id=mid, item=f'{tp} {r} {h}',
                                     curated=a, source=b, rel_diff=abs(a - b) / abs(b)))
    return rows


def stroke_volume(sheet):
    hit = sheet[sheet.apply(lambda r: r.astype(str).str.strip().str.lower().eq('total').any(), axis=1)]
    return float(hit.iloc[0].dropna().values[-1]) if len(hit) else 0.0   # 'No stroke' sheet -> 0


def check_infarct(d):
    rows = []
    for name, sheet in pd.read_excel(STROKE_XLSX, sheet_name=None, header=None).items():
        n = name.strip()
        mid = f'DBSI_M0{n.split()[-1]}' if 'DBSI' in n else f'M{int(n[2:]):02d}'
        a, b = d['bl'].loc[mid, 'stroke_size'], stroke_volume(sheet)
        rows.append(dict(check='infarct volume', mouse_id=mid, item='stroke_size (mm3)',
                         curated=a, source=b, rel_diff=abs(a - b) / max(abs(b), 1e-12) if b else abs(a)))
    return rows


def main():
    d = load_data()
    df = pd.DataFrame(check_so_power(d) + check_infarct(d))
    df['ok'] = df.rel_diff <= REL_TOL
    config.STATS.mkdir(parents=True, exist_ok=True)
    df.to_csv(config.STATS / 'input_checks.csv', index=False)
    for c, g in df.groupby('check'):
        print(f'{c:15s}: {g.ok.sum()}/{len(g)} match (max rel. diff {g.rel_diff.max():.1e})')
    if not df.ok.all():
        print(df[~df.ok].to_string(index=False))
        sys.exit(1)


if __name__ == '__main__':
    main()
