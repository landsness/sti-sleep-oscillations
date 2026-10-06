"""infarct_depth.py -- infarct depth relative to cortical thickness (Suppl Fig 1, R1 major).

Input: data/depth_traces/<mouse>_RoiSet.zip -- ImageJ tracings (J. Lee) on the
coronal cresyl-violet section of maximal infarct area, four ROIs per mouse in
this order: infarct (polygon), pia (polyline), white-matter / gray-white
boundary 'wm' (polyline), midline (polyline). Ms1 had no delineable infarct.

Normalised cortical coordinates for a point p inside the cortex:
  depth  v(p) = d(p, pia) / (d(p, pia) + d(p, wm))        0 = pia, 1 = gray/white
  medio-lateral u(p) = integral_0^s ds' / T(s')             in cortical-thickness units
     where s is the arc length along the pia from the pia end nearest the
     midline to the pia point nearest p, and T(s') is the local cortical
     thickness (distance from pia(s') to the wm line).
Infarct depth fraction = max v over the infarct outline (densified to 0.5 px;
the maximum of v over a polygon lies on its boundary), <= 0.8 = 'spares deep cortex'.

Outputs: results/tables/infarct_depth_fractions.csv  (mouse, mouse_id, STI, depth_fraction)
         results/tables/infarct_incidence_maps.npz  (u/v grid, per-group % incidence)
"""
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import roifile
from shapely.geometry import LineString, Point, Polygon
from matplotlib.path import Path as MPath

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import config  # noqa: E402
from common import load_data  # noqa: E402

ROI_ORDER = ['infarct', 'pia', 'wm', 'midline']
STEP_PX = 1.5                                  # raster step inside the infarct (image px)
U_EDGES = np.arange(0.5, 3.0001, 0.04)        # medio-lateral bins (thickness units)
V_EDGES = np.arange(0.0, 1.15001, 0.0125)      # depth bins (pia=0, gray/white=1)
DEEP_SPARING_THRESHOLD = 0.8


def mouse_to_id(name):
    """'Ms7' -> 'M07', 'DBSI_Ms2' -> 'DBSI_M02'."""
    if name.startswith('DBSI_Ms'):
        return f'DBSI_M{int(name[7:]):02d}'
    return f'M{int(name[2:]):02d}'


def read_traces(zpath):
    rois = roifile.roiread(str(zpath))
    return dict(zip(ROI_ORDER, [r.coordinates() for r in rois]))   # Ms23 ROIs are unnamed: use order


def cortical_coords(pts, pia, wm, mid):
    pia_l, wm_l, mid_l = LineString(pia), LineString(wm), LineString(mid)
    # orient the pia so that s = 0 at the end nearest the midline
    if mid_l.distance(Point(pia[0])) > mid_l.distance(Point(pia[-1])):
        pia_l = LineString(pia[::-1])
    s_grid = np.linspace(0, pia_l.length, 400)
    T = np.array([wm_l.distance(pia_l.interpolate(s)) for s in s_grid])
    U_cum = np.concatenate([[0], np.cumsum(np.diff(s_grid) / (0.5 * (T[1:] + T[:-1])))])
    u, v = [], []
    for p in pts:
        P = Point(p)
        dp, dw = pia_l.distance(P), wm_l.distance(P)
        v.append(dp / (dp + dw))
        u.append(np.interp(pia_l.project(P), s_grid, U_cum))
    return np.array(u), np.array(v)


def analyse():
    data = load_data()
    sti = data['bl']['secondary_thalamic_injury']
    rows, occ = [], {}
    for z in sorted(config.DEPTH_TRACES.glob('*_RoiSet.zip')):
        name = z.name.replace('_RoiSet.zip', '')
        R = read_traces(z)
        inf = Polygon(R['infarct']).buffer(0)
        x0, y0, x1, y1 = inf.bounds
        xs, ys = np.meshgrid(np.arange(x0, x1 + STEP_PX, STEP_PX), np.arange(y0, y1 + STEP_PX, STEP_PX))
        grid = np.c_[xs.ravel(), ys.ravel()]
        pts = grid[MPath(np.asarray(inf.exterior.coords)).contains_points(grid)]
        u, v = cortical_coords(pts, R['pia'], R['wm'], R['midline'])
        ring = inf.exterior
        edge = np.array([ring.interpolate(t).coords[0] for t in np.arange(0, ring.length, 0.5)])
        _, v_edge = cortical_coords(edge, R['pia'], R['wm'], R['midline'])
        H, _, _ = np.histogram2d(v, u, bins=[V_EDGES, U_EDGES])
        occ[name] = H > 0
        mid = mouse_to_id(name)
        rows.append(dict(mouse=name, mouse_id=mid, STI=int(sti[mid]), depth_fraction=float(v_edge.max())))
    df = pd.DataFrame(rows)

    maps = {}
    for g, label in ((1, 'STI+'), (0, 'STI-')):
        ids = sti.index[sti == g]                         # full group (mice without a traced infarct count as 0)
        stack = [occ[n] for n in occ if mouse_to_id(n) in ids]
        maps[label] = 100 * np.sum(stack, 0) / len(ids)
        maps[label + '_N'] = len(ids)
    return df, maps


def main():
    config.TABLES.mkdir(parents=True, exist_ok=True)
    df, maps = analyse()
    df.to_csv(config.TABLES / 'infarct_depth_fractions.csv', index=False)
    np.savez(config.TABLES / 'infarct_incidence_maps.npz', u_edges=U_EDGES, v_edges=V_EDGES,
             sti_pos=maps['STI+'], sti_neg=maps['STI-'], n_pos=maps['STI+_N'], n_neg=maps['STI-_N'])
    d = df.depth_fraction
    print(df.round(3).to_string(index=False))
    print(f'\n{len(d)} mice with a delineable infarct (of 25; Ms1 had none)')
    print(f'spared deep cortex (depth fraction <= {DEEP_SPARING_THRESHOLD}): {(d <= DEEP_SPARING_THRESHOLD).sum()} of {len(d)}')
    print(f'median depth fraction: {d.median():.2f}  (STI+ {df[df.STI == 1].depth_fraction.median():.2f}, '
          f'STI- {df[df.STI == 0].depth_fraction.median():.2f})')


if __name__ == '__main__':
    main()
