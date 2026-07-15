"""Supplemental Figure S3: Acute Laterality Index — 6-panel (2x3)

Row 0 (acute behavior):
  A: LI ROI1 vs acute behavior
  B: LI ROI2 vs acute behavior

Row 1 (week 1 behavior):
  C: LI ROI1 vs week 1 behavior
  D: LI ROI2 vs week 1 behavior

Row 2 (infarct volume — injury confound):
  E: LI ROI1 vs infarct volume
  F: LI ROI2 vs infarct volume

n=25, Spearman throughout.
"""
import sys
from pathlib import Path
import matplotlib.patches as mpatches
import matplotlib.pyplot as plt
from scipy.stats import spearmanr

FIGURES_DIR = Path(__file__).parent
sys.path.insert(0, str(FIGURES_DIR))

from figure_style import (
    COLORS, FONT_SIZES, add_panel_label, format_stats,
    scatter_with_regression, save_figure, set_publication_style, load_data,
)

DATA_FILE = FIGURES_DIR.parent / 'WFCI AI Summary.xlsx'
OUTPUT_DIR = FIGURES_DIR / 'final'


def main():
    set_publication_style()
    data = load_data(DATA_FILE)
    bl, ac, wk = data['bl'], data['ac'], data['wk']

    colors_all = [COLORS['sti_pos'] if bl.loc[m, 'secondary_thalamic_injury'] == 1
                  else COLORS['sti_neg'] for m in bl.index]
    patch_pos = mpatches.Patch(color=COLORS['sti_pos'], label='STI+')
    patch_neg = mpatches.Patch(color=COLORS['sti_neg'], label='STI-')

    li_r1 = ac['LI_roi1'].values
    li_r2 = ac['LI_roi2'].values
    y_ac  = ac['behavior_score'].values
    y_wk  = wk.loc[bl.index, 'behavior_score'].values
    x_inf = bl['stroke_size'].values

    fig, axes = plt.subplots(3, 2, figsize=(10, 13))

    # Per-panel spec: (ax, x, y, xlabel, ylabel, stats_loc, legend_loc, zero_line, panel)
    # Stats and legend locations chosen to avoid data overlap for each correlation direction.
    # All behavior panels have negative correlation (upper-right is sparse) — stats upper right.
    # Legend (Panel A only) goes lower right to avoid the stats box.
    # Infarct panels also negative — stats upper right is fine (no legend).
    panels = [
        (axes[0, 0], li_r1, y_ac,  'Acute Laterality Index (Lesion)',
         'Acute Forelimb Asymmetry Score (%)',  'upper right', 'upper left', True,  'A'),
        (axes[0, 1], li_r2, y_ac,  'Acute Laterality Index (Perilesional)',
         'Acute Forelimb Asymmetry Score (%)',  'upper right', None,          True,  'B'),
        (axes[1, 0], li_r1, y_wk,  'Acute Laterality Index (Lesion)',
         'Week 1 Forelimb Asymmetry Score (%)', 'upper right', None,          True,  'C'),
        (axes[1, 1], li_r2, y_wk,  'Acute Laterality Index (Perilesional)',
         'Week 1 Forelimb Asymmetry Score (%)', 'upper right', None,          True,  'D'),
        (axes[2, 0], li_r1, x_inf, 'Acute Laterality Index (Lesion)',
         'Infarct Volume (mm³)',                'upper right', None,          False, 'E'),
        (axes[2, 1], li_r2, x_inf, 'Acute Laterality Index (Perilesional)',
         'Infarct Volume (mm³)',                'upper right', None,          False, 'F'),
    ]

    for ax, x, y, xlabel, ylabel, stats_loc, legend_loc, zero_line, panel in panels:
        rho, p = spearmanr(x, y)
        scatter_with_regression(
            ax, x, y, colors_all,
            xlabel=xlabel, ylabel=ylabel,
            stats_text=format_stats(rho=rho, p=p),
            loc=stats_loc,
        )
        if zero_line:
            ax.axhline(0, color='#666666', lw=0.8, ls='--', alpha=0.6, zorder=1)
        if legend_loc:
            ax.legend(handles=[patch_pos, patch_neg],
                      fontsize=FONT_SIZES['legend'], framealpha=0.8, loc=legend_loc)
        add_panel_label(ax, panel)

    plt.tight_layout(pad=0.5)
    save_figure(fig, 'SupplFig_S3_AcuteLI', OUTPUT_DIR)
    plt.close()
    print('SupplFig S3 done.')


if __name__ == '__main__':
    main()
