#!/usr/bin/env python3
"""
fig2b_lesional_traces.py  (was plot_roi_traces.py)
Lesional-ROI slow-oscillation traces (baseline / acute / week 1), overlaying
uncorrected vs hemodynamically corrected dF/F. Answers the reviewer request to
"show hemodynamically corrected fluorescence traces."

Input : data/derived/ROI_Traces_export.mat  (from upstream/matlab/Extract_ROI_Traces.m)
Output: results/panels/Fig2B_lesional_traces.{png,pdf,svg}   (stand-alone review panel: baseline + acute)
        python fig2b_lesional_traces.py --all-sessions  adds week 1
Figure 2 itself calls draw(axs, ...) to render these traces directly in its own axes
(rev. 2026-10-04), so fonts match the rest of the figure at print size.

Representative mouse chosen by SO signal-to-noise (SO 0.1-1 Hz / 1-5 Hz power)
on the baseline lesional corrected trace; set MOUSE manually to override.
"""
import sys
from pathlib import Path
import numpy as np
from scipy.io import loadmat
from scipy import signal
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
import config  # noqa: E402
INPUT = config.ROI_TRACES
MOUSE = 8               # locked representative (highest SO SNR among near-median mice 8/9/14/24); None -> auto-pick
ROI   = "lesional"
WIN   = 60.0            # seconds to plot
C_CORR, C_UNC = "#000000", "#A6A6A6"   # corrected black, uncorrected light gray (orange is reserved for "acute" in 2C)
SESSIONS = [('bsl', 'Baseline'), ('acute', 'Acute (24 h)'), ('wk1', 'Week 1')]


def load():
    d = loadmat(INPUT, squeeze_me=True, struct_as_record=False)
    D = {}
    for e in d["EXPORT"]:
        D[(int(e.mouse), str(e.session), str(e.roi), str(e.corr))] = (np.asarray(e.trace, float).ravel(), float(e.fs))
    return D


def snr(tr, fs):
    tr = tr - np.nanmean(tr)
    f, P = signal.welch(tr, fs=fs, nperseg=min(256, len(tr)))
    return P[(f >= 0.1) & (f <= 1)].mean() / P[(f > 1) & (f <= 5)].mean()


def pick_mouse(D):
    if MOUSE is not None:
        return MOUSE
    mice = sorted({k[0] for k in D})
    scores = {m: snr(*D[(m, 'bsl', ROI, 'corr')]) for m in mice if (m, 'bsl', ROI, 'corr') in D}
    m = max(scores, key=scores.get)
    print("SNR:", {k: round(v, 1) for k, v in scores.items()}, "-> Ms%d" % m)
    return m


def draw(axs, sessions=SESSIONS[:2], fs=None, lw=(0.8, 1.0)):
    """Draw stacked traces into axs (one per session). fs: dict(label, tick, legend, title)."""
    fs = fs or dict(label=10, tick=9, legend=9, title=11)
    D = load(); m = pick_mouse(D)
    allv = []
    for s, _ in sessions:
        for c in ('uncorr', 'corr'):
            if (m, s, ROI, c) in D:
                tr, f = D[(m, s, ROI, c)]; n = int(WIN * f); allv.append(tr[:n] - np.mean(tr[:n]))
    ylim = np.nanpercentile(np.abs(np.concatenate(allv)), 99.5) * 1.15
    for ai, (s, lab) in enumerate(sessions):
        ax = axs[ai]
        for c, col, w, name in [('uncorr', C_UNC, lw[0], 'Uncorrected'), ('corr', C_CORR, lw[1], 'Corrected')]:
            if (m, s, ROI, c) not in D:
                continue
            tr, f = D[(m, s, ROI, c)]; n = int(WIN * f); t = np.arange(n) / f
            ax.plot(t, tr[:n] - np.mean(tr[:n]), color=col, lw=w, label=name)
        ax.set_ylim(-ylim, ylim); ax.set_xlim(0, WIN)
        ax.set_ylabel(r'$\Delta$F/F', fontsize=fs['label'], labelpad=1)
        ax.tick_params(labelsize=fs['tick'], length=2.5, pad=1.5)
        ax.text(0.01, 0.97, lab, transform=ax.transAxes, fontsize=fs['title'], fontweight='bold', va='top')
        ax.axhline(0, color='#bbb', lw=0.5, zorder=0)
        for sp in ('top', 'right'):
            ax.spines[sp].set_visible(False)
        if ai < len(sessions) - 1:
            ax.tick_params(labelbottom=False)
        if ai == 0:
            ax.legend(frameon=False, fontsize=fs['legend'], loc='lower right', bbox_to_anchor=(1.0, 1.0), ncol=2,
                      handlelength=1.4, columnspacing=1.0, borderaxespad=0.0)
    axs[-1].set_xlabel('Time (s)', fontsize=fs['label'], labelpad=1)
    return m


def main(all_sessions=False):
    plt.rcParams.update({"font.family": "sans-serif", "font.size": 10, "axes.linewidth": 0.8,
                         "svg.fonttype": "none", "pdf.fonttype": 42})
    sessions = SESSIONS if all_sessions else SESSIONS[:2]
    fig, axs = plt.subplots(len(sessions), 1, figsize=(9, 7 * len(sessions) / 3), sharex=True, sharey=True)
    m = draw(axs, sessions, lw=(1.0, 1.3))
    if all_sessions:   # stand-alone review version keeps its own sub-panel letters
        for ai, ax in enumerate(axs):
            ax.text(-0.02, 1.02, "ABC"[ai], transform=ax.transAxes, fontsize=14, fontweight='bold', va='bottom', ha='right')
        fig.suptitle(f'Lesional ROI slow-oscillation traces (Ms{m}) — corrected vs uncorrected ΔF/F', fontsize=12, y=0.98)
    plt.tight_layout(rect=[0, 0, 1, 0.96] if all_sessions else None)
    config.PANELS.mkdir(parents=True, exist_ok=True)
    stem = 'Fig2B_lesional_traces' + ('_all_sessions' if all_sessions else '')
    for ext in ('png', 'pdf', 'svg'):
        fig.savefig(config.PANELS / f'{stem}.{ext}', dpi=300, bbox_inches='tight')
    print("saved", config.PANELS / f"{stem}.png", "(Ms%d)" % m)


if __name__ == '__main__':
    main('--all-sessions' in sys.argv)
