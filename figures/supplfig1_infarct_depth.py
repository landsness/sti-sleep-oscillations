"""supplfig1_infarct_depth.py -- Supplementary Figure 1: coronal lesion-incidence maps by STI group.

Group maps in a normalised cortical frame (depth: pia = 0, gray/white = 1;
medio-lateral: distance from midline in cortical-thickness units), colour =
% of animals in the group with infarct at that location. Computed by
analysis/infarct_depth.py from J. Lee's ImageJ tracings (data/depth_traces/).

Output: results/figures/{png,pdf,svg}/SupplFig1_InfarctDepth.*
"""
import sys
from pathlib import Path

import numpy as np
import matplotlib.pyplot as plt

FIGURES_DIR = Path(__file__).resolve().parent
sys.path[:0] = [str(FIGURES_DIR), str(FIGURES_DIR.parent), str(FIGURES_DIR.parent / 'analysis')]
import config  # noqa: E402
from figure_style import set_publication_style, save_figure  # noqa: E402
import infarct_depth  # noqa: E402


def main():
    df, maps = infarct_depth.analyse()
    set_publication_style()
    ue, ve = infarct_depth.U_EDGES, infarct_depth.V_EDGES
    fig, axes = plt.subplots(1, 2, figsize=(10, 3.4), sharey=True)
    for ax, label, g in zip(axes, ('STI+', 'STI-'), (1, 0)):
        med = df[df.STI == g].depth_fraction.median()
        im = ax.pcolormesh(ue, ve, maps[label], cmap='magma', vmin=0, vmax=100, shading='flat')
        ax.axhline(1.0, color='white', ls='--', lw=1)
        ax.axhline(0.0, color='white', lw=0.8)
        ax.set_ylim(ve[-1], 0)
        ax.set_aspect('equal')
        ax.set_title(f"{label.replace('-', '−')}  (N={maps[label + '_N']}, median depth {med:.2f})", fontsize=9)
        ax.set_xlabel('Distance from midline (cortical-thickness units)', fontsize=8)
        if g == 1:
            ax.set_ylabel('Cortical depth\n(pia=0, WM=1)', fontsize=8)
            ax.text(ue[-1] - 0.05, 0.05, 'pia', color='white', ha='right', va='top', fontsize=7)
            ax.text(ue[-1] - 0.05, 0.97, 'gray/white', color='white', ha='right', va='bottom', fontsize=7)
    cb = fig.colorbar(im, ax=axes, fraction=0.025, pad=0.02)
    cb.set_label('% of mice with infarct', fontsize=8)
    save_figure(fig, 'SupplFig1_InfarctDepth', config.FIG_OUT)
    plt.close(fig)


if __name__ == '__main__':
    main()
