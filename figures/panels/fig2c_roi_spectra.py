#!/usr/bin/env python3
"""
fig2c_roi_spectra.py  (was plot_roi_spectra.py)
Group power-spectra figure (3 panels: lesional / perilesional / contralateral),
group mean +/- SEM across mice, log-log, 0.1-1 Hz SO band shaded.

Input : data/derived/ROI_Spectra_export.mat  (from upstream/matlab/Extract_ROI_Spectra.m)
Output: results/panels/Fig2C_roi_spectra.{png,pdf,svg}   (stand-alone review panel: baseline + acute)
        python fig2c_roi_spectra.py --all-sessions  adds week 1
Figure 2 itself calls draw(axes, ...) to render the spectra directly in its own axes
(rev. 2026-10-04), so fonts match the rest of the figure at print size.

Normalization: each mouse's spectra are divided by that mouse's baseline
contralateral broadband geometric-mean power (0.02-10 Hz), removing per-animal
gain so the group mean is meaningful. Curves are geometric mean; band = +/-SEM
in log10 space.

Cohort (rev. 2026-10-04): baseline and acute use all 25 animals; week 1 omits
the two pilot animals in config.WEEK1_SO_EXCLUDED_ROI_IDS (n = 23), matching
every other SO analysis. n is counted from the data, not hard-coded.
The ~8 Hz heart-rate peak is explained in the figure legend (no in-panel callout).
"""
import sys
from pathlib import Path
import numpy as np
from scipy.io import loadmat
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
import config  # noqa: E402
INPUT = config.ROI_SPECTRA
FMIN, FMAX = 0.02, 10.0            # plotted/normalization frequency range
SO_LO, SO_HI = 0.1, 1.0            # slow-oscillation band
ROIS = ["lesional", "perilesional", "contra"]
TITLES = {"lesional": "Lesional", "perilesional": "Perilesional",
          "contra": "Contralateral homolog"}
LABELS = {"bsl": "Baseline", "acute": "Acute (24 h)", "wk1": "Week 1"}
COLORS = {"bsl": "#000000", "acute": "#D55E00", "wk1": "#0072B2"}   # Okabe-Ito


def load():
    d = loadmat(INPUT, squeeze_me=True, struct_as_record=False)
    hz = np.asarray(d["hz"]).ravel()
    bym = {}
    excl = set(config.WEEK1_SO_EXCLUDED_ROI_IDS)
    for e in d["EXPORT"]:
        m = str(e.mouse).strip()           # '7', '2DBSI', ... (older exports stored ints)
        if str(e.session) == "wk1" and m in excl:
            continue                       # week-1 SO excluded for the MRI pilot animals
        bym.setdefault(m, {})[(str(e.roi), str(e.session))] = np.asarray(e.psd, float).ravel()
    mice = sorted(bym, key=lambda s: (not s.isdigit(), int(s) if s.isdigit() else s))
    fmask = (hz >= FMIN) & (hz <= FMAX)
    NF = {}
    for m in mice:
        for key in [("contra", "bsl"), ("lesional", "bsl"), ("perilesional", "bsl")]:
            if key in bym[m]:
                v = bym[m][key][fmask]; v = v[v > 0]
                if v.size:
                    NF[m] = np.exp(np.mean(np.log(v))); break

    def stack(roi, sess):
        return np.asarray([bym[m][(roi, sess)][fmask] / NF[m]
                           for m in mice if (roi, sess) in bym[m] and NF.get(m)])
    return hz[fmask], stack, mice


def draw(axes, sess=("bsl", "acute"), fs=None, lw=1.8, legend_ax=0):
    """Draw the three ROI spectra into axes (len 3, shared y). fs: dict(label, tick, legend, title, note)."""
    fs = fs or dict(label=10, tick=9, legend=9, title=11, note=8)
    f, stack, mice = load()
    print(f"fig2c: {len(mice)} mice in export")
    ns = {}
    for ai, roi in enumerate(ROIS):
        ax = axes[ai]
        ax.axvspan(SO_LO, SO_HI, color="#F5D400", alpha=0.13, lw=0, zorder=0)
        for s in sess:
            A = stack(roi, s)
            if A.size == 0:
                continue
            logA = np.log10(A)
            mu = logA.mean(0); sem = logA.std(0, ddof=1) / np.sqrt(A.shape[0])
            ax.fill_between(f, 10**(mu - sem), 10**(mu + sem), color=COLORS[s], alpha=0.18, lw=0)
            ax.plot(f, 10**mu, color=COLORS[s], lw=lw, label=f"{LABELS[s]} (n = {A.shape[0]})")
            ns[s] = A.shape[0]
        ax.set_xscale("log"); ax.set_yscale("log")
        ax.set_xlim(FMIN, FMAX)
        ax.set_xlabel("Frequency (Hz)", fontsize=fs['label'], labelpad=1)
        ax.set_title(TITLES[roi], fontsize=fs['title'], pad=3)
        ax.tick_params(labelsize=fs['tick'], length=2.5, pad=1.5)
        for sp in ('top', 'right'):
            ax.spines[sp].set_visible(False)
        ax.text(np.sqrt(SO_LO * SO_HI), 0.99, "SO band", transform=ax.get_xaxis_transform(),
                ha="center", va="top", fontsize=fs['note'], color="#8a6d00")
        if ai == 0:
            ax.set_ylabel("Normalized power (a.u.)", fontsize=fs['label'], labelpad=1)
        if ai == legend_ax:
            ax.legend(frameon=False, fontsize=fs['legend'], loc="lower left", handlelength=1.4, borderaxespad=0.2)
    return ns


def main(all_sessions=False):
    plt.rcParams.update({"font.family": "sans-serif", "font.size": 10, "axes.linewidth": 0.8,
                         "svg.fonttype": "none", "pdf.fonttype": 42, "ps.fonttype": 42})
    sess = ("bsl", "acute", "wk1") if all_sessions else ("bsl", "acute")
    fig, axes = plt.subplots(1, 3, figsize=(13, 4.6), sharey=True)
    ns = draw(axes, sess)
    if all_sessions:   # stand-alone review version keeps its own sub-panel letters and title
        for ai, ax in enumerate(axes):
            ax.text(-0.02, 1.04, "ABC"[ai], transform=ax.transAxes, fontsize=15, fontweight="bold", va="bottom", ha="right")
        fig.suptitle(f"ROI power spectra (group mean ± SEM; n = {ns['bsl']}, week 1 n = {ns['wk1']}) "
                     "— hemodynamically corrected calcium", fontsize=12, y=0.99)
    plt.tight_layout(rect=[0, 0, 1, 0.96] if all_sessions else None)
    config.PANELS.mkdir(parents=True, exist_ok=True)
    stem = config.PANELS / ("Fig2C_roi_spectra" + ("_all_sessions" if all_sessions else ""))
    for ext in ("png", "pdf", "svg"):
        fig.savefig(f"{stem}.{ext}", dpi=300, bbox_inches="tight")
    print("saved:", f"{stem}.{{png,pdf,svg}}")


if __name__ == "__main__":
    main("--all-sessions" in sys.argv)
