"""fig2a_lesion_topo.py -- Figure 2A: bilateral group map of acute SO power (% baseline).

Rev. 2026-10-04: redesigned from the lesion-centred, ipsilesional-only crop to a
bilateral map so the panel also serves as the field-of-view / ROI schematic
requested by Reviewer 2 (minor 1). Output stem unchanged (Fig2A_lesion_topo).

Inputs
  data/rawmaps/Mouse_<ID>_SO_RawMaps.mat   so_map_bsl / so_map_acute / so_map_wk1
      128x128 SO-band (0.1-1 Hz) mean PSD in the bregma-centred atlas frame,
      written by upstream/matlab/SO_PowerAnalysis_Arnav_v2.m (hemodynamically
      corrected calcium, contrast 4 of the fc1 Power.mat).
  data/ROI_Summary_V2.xlsx                  ROI1/ROI2 centroids and pixel counts
      (SO_ROI_Quantification_mirrored_v2.m). The drawn ROI masks were not saved,
      so ROIs are overlaid as equal-area circles at the cohort-mean centroid;
      contralateral ROIs are their mirror images across the midline (x -> 129 - x).
  data/Fig2A_display_outline.csv            closed outline of the displayed cortical field
      (x, y vertices, 0-indexed pixel coords). Traced by E.L. on the acute group map
      (2026-10-04), then made mirror-symmetric about the midline (mean of each side
      and the mirrored other side, in polar coordinates) and smoothed (Fourier series
      truncated at 10 harmonics). Display only: it hides window-edge pixels and does
      not affect any ROI value or statistic.

Method
  1. Per mouse, pct = 100 * so_map_<tp> / so_map_bsl inside the brain mask
     (so_map_bsl > 5% of its 99th percentile; the max was used
     before 2026-10-04 and wrongly dropped valid cortex in Mouse 2, whose map has one
     very bright spot). No re-centring: maps stay in the bregma-
     registered atlas frame.
  2. Animals with a right-hemisphere lesion (ROI1 centroid x > 64.5; 2DBSI) are
     mirrored left-right so the ipsilesional hemisphere is always on the left.
  3. Average across mice (pixels valid in >= 50% of mice). Inside the display outline,
     the few edge pixels below that coverage are filled from the nearest valid pixel
     (display only); everything outside the outline is blank.
  4. Overlay lesional (solid) and perilesional (dashed) ROIs and their mirrors,
     the 1-mm photothrombotic spot (dotted) and a 1-mm scale bar.
Cohort: acute n = 25 (all animals); week 1 (--week1) n = 23, omitting
config.WEEK1_SO_EXCLUDED_ROI_IDS.

Outputs: results/panels/Fig2A_lesion_topo.{png,pdf,svg}
         (--week1 also writes Fig2A_lesion_topo_acute_wk1.* with a week-1 panel)
Sanity check printed: group value at the four ROI centres.
"""
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.io import loadmat
from scipy import ndimage
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patheffects as pe
from matplotlib.lines import Line2D
from matplotlib.path import Path as MplPath
from matplotlib.patches import PathPatch

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
import config  # noqa: E402

ALL_MICE = [str(i) for i in range(1, 25) if i != 18] + ['2DBSI', '4DBSI']   # n = 25
MID = 63.5                      # midline, 0-indexed (64.5 in MATLAB coordinates)
COVERAGE = 0.50
CLIM = (0, 100)                 # map spans 14-97 %; nothing exceeds baseline
CMAP = 'magma'                  # perceptually uniform, CVD- and greyscale-safe (chosen 2026-10-04)
DISPLAY_OUTLINE = config.DATA / 'Fig2A_display_outline.csv'


def load_maps():
    roi = pd.read_excel(config.ROI_SUMMARY)
    roi['Mouse'] = roi['Mouse'].astype(str)
    roi = roi.set_index('Mouse')
    wk1_mice = [m for m in ALL_MICE if m not in set(config.WEEK1_SO_EXCLUDED_ROI_IDS)]
    stacks = {'acute': [], 'wk1': []}
    cents = []
    for m in ALL_MICE:
        d = loadmat(config.RAWMAPS / f'Mouse_{m}_SO_RawMaps.mat')
        bsl = d['so_map_bsl'].astype(float)
        brain = np.isfinite(bsl) & (bsl > 0.05 * np.nanpercentile(bsl, 99))   # robust to single bright pixels (Mouse 2)
        r = roi.loc[m]
        x1, y1 = r.ROI1_CentroidX - 1, r.ROI1_CentroidY - 1      # 0-indexed
        x2, y2 = r.ROI2_CentroidX - 1, r.ROI2_CentroidY - 1
        flip = r.ROI1_CentroidX > 64.5                           # right-hemisphere lesion
        if flip:
            x1, x2 = 127 - x1, 127 - x2
        cents.append((x1, y1, x2, y2, r.ROI1_Pixels, r.ROI2_Pixels))
        for tp in stacks:
            if tp == 'wk1' and m not in wk1_mice:
                continue
            pct = 100 * d[f'so_map_{tp}'].astype(float) / np.where(bsl > 0, bsl, np.nan)
            pct[~brain] = np.nan
            stacks[tp].append(pct[:, ::-1] if flip else pct)
    c = np.array(cents, float)
    geom = dict(x1=c[:, 0].mean(), y1=c[:, 1].mean(), x2=c[:, 2].mean(), y2=c[:, 3].mean(),
                r1=np.sqrt(c[:, 4].mean() / np.pi), r2=np.sqrt(c[:, 5].mean() / np.pi),
                rL=config.PT_SPOT_RADIUS_MM / config.MM_PER_PX)
    outline = MplPath(np.loadtxt(DISPLAY_OUTLINE, delimiter=','), closed=True)
    yy, xx = np.mgrid[0:128, 0:128]
    inside = outline.contains_points(np.c_[xx.ravel(), yy.ravel()]).reshape(128, 128)
    groups, ns = {}, {}
    for tp, st in stacks.items():
        A = np.stack(st, 2)
        G = np.nanmean(A, 2)
        G[np.isfinite(A).mean(2) < COVERAGE] = np.nan
        # fill low-coverage pixels inside the outline from the nearest valid pixel
        idx = ndimage.distance_transform_edt(~np.isfinite(G), return_distances=False, return_indices=True)
        G = G[tuple(idx)]
        G[~ndimage.binary_dilation(inside, iterations=2)] = np.nan   # 2-px margin; the vector outline clips the edge
        groups[tp], ns[tp] = G, A.shape[2]
    return groups, ns, geom, outline


def _circle(ax, x, y, r, ls, lw=1.4):
    t = np.linspace(0, 2 * np.pi, 200)
    ax.plot(x + r * np.cos(t), y + r * np.sin(t), color='black', ls=ls, lw=lw,
            path_effects=[pe.withStroke(linewidth=lw + 1.6, foreground='white')])


FS = dict(title=10, hemi=8, note=6.5, scale=7, legend=6.5)     # stand-alone sizes; Figure 2 passes its own


def draw(ax, G, geom, outline, title=None, legend=False, scalebar=True, hemi_labels=True,
         fs=None, ring_lw=1.4):
    fs = {**FS, **(fs or {})}
    cmap = plt.get_cmap(CMAP).copy()
    cmap.set_bad('white')
    im = ax.imshow(G, cmap=cmap, vmin=CLIM[0], vmax=CLIM[1], interpolation='nearest')
    im.set_clip_path(PathPatch(outline, transform=ax.transData, facecolor='none'))
    ax.add_patch(PathPatch(outline, facecolor='none', edgecolor='0.55', lw=0.6))
    for x, y, r, ls in ((geom['x1'], geom['y1'], geom['r1'], '-'),
                        (geom['x2'], geom['y2'], geom['r2'], '--')):
        _circle(ax, x, y, r, ls, lw=ring_lw)
        _circle(ax, 127 - x, y, r, ls, lw=ring_lw)                            # contralateral mirror
    _circle(ax, geom['x1'], geom['y1'], geom['rL'], ':', lw=ring_lw * 0.85)
    v = outline.vertices
    xmin, xmax, ymin, ymax = v[:, 0].min(), v[:, 0].max(), v[:, 1].min(), v[:, 1].max()
    ax.set_xlim(xmin - 3, xmax + 3)
    ax.set_ylim(ymax + 6, ymin - 9)
    ax.axis('off')
    if title:
        ax.set_title(title, fontsize=fs['title'], fontweight='bold', pad=2)
    if hemi_labels:
        for x, txt in ((MID - 30, 'Ipsilesional'), (MID + 30, 'Contralesional')):
            ax.text(x, ymin - 4, txt, ha='center', va='bottom', fontsize=fs['hemi'], fontweight='bold')
        ax.annotate('', xy=(xmin - 1, ymin + 1), xytext=(xmin - 1, ymin + 12),
                    arrowprops=dict(arrowstyle='-|>', lw=0.9, color='0.25', mutation_scale=7))
        ax.text(xmin + 1.5, ymin + 3, 'Anterior', fontsize=fs['note'], color='0.25', ha='left', va='center')
    if scalebar:
        L = 1.0 / config.MM_PER_PX
        xb, yb = xmin + L, ymax + 3.5
        ax.plot([xb - L, xb], [yb, yb], color='black', lw=2.2, solid_capstyle='butt')
        ax.text(xb - L / 2, yb - 1.2, '1 mm', ha='center', va='bottom', fontsize=fs['scale'])
    if legend:
        h = [Line2D([], [], color='k', ls='-', label='Lesional ROI'),
             Line2D([], [], color='k', ls='--', label='Perilesional ROI'),
             Line2D([], [], color='k', ls=':', label='1-mm PT spot')]
        ax.legend(handles=h, loc='upper center', bbox_to_anchor=(0.5, 0.0), ncol=3, fontsize=fs['legend'],
                  frameon=False, handlelength=2.0, columnspacing=0.9, handletextpad=0.5)
    return im


def _save(fig, stem):
    config.PANELS.mkdir(parents=True, exist_ok=True)
    for ext in ('png', 'pdf', 'svg'):
        fig.savefig(config.PANELS / f'{stem}.{ext}', dpi=300, bbox_inches='tight')
    print('saved', config.PANELS / f'{stem}.png')


def main(week1=False):
    plt.rcParams.update({'font.family': 'sans-serif', 'svg.fonttype': 'none', 'pdf.fonttype': 42})
    groups, ns, geom, outline = load_maps()
    for tp, G in groups.items():
        vals = [G[int(round(geom['y1'])), int(round(geom['x1']))],
                G[int(round(geom['y2'])), int(round(geom['x2']))],
                G[int(round(geom['y1'])), int(round(127 - geom['x1']))],
                G[int(round(geom['y2'])), int(round(127 - geom['x2']))]]
        print(f'{tp:5s} (n={ns[tp]}): lesional {vals[0]:.1f}% | perilesional {vals[1]:.1f}% | '
              f'contra lesional {vals[2]:.1f}% | contra perilesional {vals[3]:.1f}% of baseline')

    fig, ax = plt.subplots(figsize=(3.9, 3.3))
    im = draw(ax, groups['acute'], geom, outline, title='Acute (24 h)', legend=True)
    cb = fig.colorbar(im, ax=ax, fraction=0.04, pad=0.02, extend='neither')
    cb.set_label('SO power (% baseline)', fontsize=8)
    cb.ax.tick_params(labelsize=7)
    _save(fig, 'Fig2A_lesion_topo')
    plt.close(fig)

    if week1:
        fig, axs = plt.subplots(1, 2, figsize=(7.6, 3.3))
        draw(axs[0], groups['acute'], geom, outline, f"Acute (24 h), n = {ns['acute']}", legend=True)
        im = draw(axs[1], groups['wk1'], geom, outline, f"Week 1, n = {ns['wk1']}")
        cb = fig.colorbar(im, ax=axs, fraction=0.03, pad=0.02, extend='neither')
        cb.set_label('SO power (% baseline)', fontsize=8)
        _save(fig, 'Fig2A_lesion_topo_acute_wk1')
        plt.close(fig)


if __name__ == '__main__':
    main(week1='--week1' in sys.argv)
