#!/usr/bin/env python3
"""
so_contrast_maps.py  — MATLAB-free version of SO_Contrast_SO_Maps.m

Extracts 0.1-1 Hz slow-oscillation (SO) power maps for ALL contrasts
(HbO, HbR, HbT, Calcium) from the existing *-Power.mat files, so you can
test whether the lesional SO "suppression" seen in calcium also appears in
the hemodynamic channels.

Reuses the same `whole_spectra_map` the calcium SO maps came from, so the
method is identical to the lab pipeline by construction. No reprocessing of
dataFluor / dataHb needed.

whole_spectra_map(:,:,freq,contrast) contrast order (LeeWrapper 301-305):
    1 = HbO   2 = HbR   3 = HbT   4 = Calcium (hemodynamically corrected)

REQUIREMENTS:  pip install h5py numpy scipy matplotlib
RUN:           edit BASE_DIR below, then:  python so_contrast_maps.py
"""

import os, glob, re
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

# ---------------- USER INPUT ----------------
# Folder that contains the per-mouse subfolders (1, 2, ... each with bsl/acute/oneweek)
BASE_DIR = r"C:\Users\landsness\Box\DBSI Directory\WFCI\WF_Calcium_RepositoryREDONE\PreProcessed_Data"
OUT_DIR  = r"C:\Users\landsness\Box\DBSI Directory\WFCI\WF_Calcium_RepositoryREDONE\SO_Analysis\Contrast_SO_Maps"
SO_LOW, SO_HIGH = 0.1, 1.0
CONTRASTS = ["HbO", "HbR", "HbT", "Calcium"]     # channels 1..4 (MATLAB) -> index 0..3
SESSION_DIRS = {"bsl": "bsl", "acute": "acute", "wk1": "oneweek"}  # logical -> folder name
RATIO_CLIM = (0, 2)

os.makedirs(OUT_DIR, exist_ok=True)


def load_power_mat(path):
    """Return (whole_spectra_map[Y,X,freq,contrast], hz) from a v7.3 or older .mat."""
    hz = None; wsm = None
    try:
        import h5py
        with h5py.File(path, "r") as f:
            hz = np.array(f["hz"]).squeeze().astype(float)
            wsm = np.array(f["whole_spectra_map"])   # h5py gives reversed dim order
    except Exception:
        from scipy.io import loadmat
        d = loadmat(path)
        hz = np.asarray(d["hz"]).squeeze().astype(float)
        wsm = np.asarray(d["whole_spectra_map"])
    # ---- reorder to [Y, X, freq, contrast] robustly, regardless of storage order ----
    L = hz.size
    shape = wsm.shape
    assert wsm.ndim == 4, f"expected 4-D whole_spectra_map, got {shape}"
    axes = list(range(4))
    freq_axis = int(np.argmin([abs(s - L) for s in shape]))            # axis whose size == len(hz)
    spatial = [a for a in axes if shape[a] == 128 and a != freq_axis][:2]
    contrast_axis = [a for a in axes if a not in spatial and a != freq_axis][0]
    wsm = np.transpose(wsm, (spatial[0], spatial[1], freq_axis, contrast_axis))
    return wsm, hz


def so_map(wsm, hz, contrast_idx):
    band = (hz > SO_LOW) & (hz < SO_HIGH)
    return np.nanmean(wsm[:, :, band, contrast_idx], axis=2)   # linear SO power


def find_power_file(session_dir):
    hits = sorted(glob.glob(os.path.join(session_dir, "*Power.mat")))
    return hits[0] if hits else None


# ---------------- MAIN ----------------
mice = sorted([d for d in os.listdir(BASE_DIR)
               if os.path.isdir(os.path.join(BASE_DIR, d)) and re.fullmatch(r"\d+|.*DBSI", d)])
print(f"Found {len(mice)} mouse folders")

for m in mice:
    SO = {c: {} for c in CONTRASTS}
    got = False
    for sess, folder in SESSION_DIRS.items():
        sd = os.path.join(BASE_DIR, m, folder)
        if not os.path.isdir(sd):
            continue
        pf = find_power_file(sd)
        if not pf:
            print(f"  [{m}/{sess}] no Power.mat"); continue
        try:
            wsm, hz = load_power_mat(pf)
        except Exception as e:
            print(f"  [{m}/{sess}] load failed: {e}"); continue
        for ci, c in enumerate(CONTRASTS):
            SO[c][sess] = so_map(wsm, hz, ci)
        got = True
        print(f"  [{m}/{sess}] ok  (freq bins={hz.size})")
    if not got:
        continue

    # save compact arrays (small — easy to share/stage back)
    np.savez_compressed(os.path.join(OUT_DIR, f"Mouse_{m}_SO_AllContrasts.npz"),
                        **{f"{c}_{s}": SO[c][s] for c in CONTRASTS for s in SO[c]})

    # ratio-map comparison figure: rows=contrast, cols= acute/bsl, wk1/bsl
    fig, ax = plt.subplots(len(CONTRASTS), 2, figsize=(6, 2.4*len(CONTRASTS)))
    for r, c in enumerate(CONTRASTS):
        for k, sess in enumerate(["acute", "wk1"]):
            a = ax[r, k]
            if "bsl" in SO[c] and sess in SO[c]:
                b = SO[c]["bsl"].copy(); b[b == 0] = np.nan
                mm = SO[c][sess].copy(); mm[mm == 0] = np.nan
                im = a.imshow(mm / b, cmap="jet", vmin=RATIO_CLIM[0], vmax=RATIO_CLIM[1])
            a.set_xticks([]); a.set_yticks([]); a.set_title(f"{c}  {sess}/bsl", fontsize=9)
    fig.suptitle(f"Mouse {m} — SO (0.1-1 Hz) ratio maps by contrast", fontsize=11)
    fig.colorbar(im, ax=ax, fraction=0.03, pad=0.02)
    fig.savefig(os.path.join(OUT_DIR, f"Mouse_{m}_SO_ContrastRatios.png"),
                dpi=150, bbox_inches="tight")
    plt.close(fig)

print(f"\nDone. Outputs (small .npz + PNGs) in:\n  {OUT_DIR}")
print("Compare the HbT/HbO/HbR ratio maps against Calcium at the lesion.")
