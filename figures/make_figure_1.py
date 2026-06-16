"""make_figure_1.py â€" Figure 1: Injury Severity Defines the Physiological Range

Three-panel strip+box plot comparing STI+ vs STI- groups for:
  A: Infarct Volume (mmÂ³)
  B: Acute Behavior Score
  C: Week 1 Behavior Score
"""

import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from scipy.stats import mannwhitneyu

# ---------------------------------------------------------------------------
# Path setup â€" allow import from figures/ regardless of working directory
# ---------------------------------------------------------------------------
FIGURES_DIR = Path(__file__).parent
sys.path.insert(0, str(FIGURES_DIR))

from figure_style import (
    COLORS, set_publication_style, strip_box_panel, save_figure, load_data,
)

# ---------------------------------------------------------------------------
# Paths
# ---------------------------------------------------------------------------
DATA_FILE = FIGURES_DIR.parent / 'WFCI AI Summary.xlsx'
OUTPUT_DIR = FIGURES_DIR / 'final'

# ---------------------------------------------------------------------------
# Load data
# ---------------------------------------------------------------------------
def load_data_fig():
    return load_data(DATA_FILE)


def prepare_data(data):
    bl, ac, wk = data['bl'], data['ac'], data['wk']
    sti_pos, sti_neg = data['sti_pos'], data['sti_neg']

    # Panel A: Infarct volume (baseline)
    stroke_pos = bl.loc[sti_pos, 'stroke_size'].values
    stroke_neg = bl.loc[sti_neg, 'stroke_size'].values

    # Panel B: Baseline behavior score (U=82, p=0.805 — null result)
    bl_beh_pos = bl.loc[sti_pos, 'behavior_score'].values
    bl_beh_neg = bl.loc[sti_neg, 'behavior_score'].values

    # Panel C: Acute (24h) behavior score (U=106, p=0.119 — trend only)
    ac_beh_pos = ac.loc[sti_pos, 'behavior_score'].values
    ac_beh_neg = ac.loc[sti_neg, 'behavior_score'].values

    return (stroke_pos, stroke_neg), (bl_beh_pos, bl_beh_neg), (ac_beh_pos, ac_beh_neg)


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------
def main():
    set_publication_style()

    data = load_data_fig()
    (stroke_pos, stroke_neg), (bl_beh_pos, bl_beh_neg), (ac_beh_pos, ac_beh_neg) = prepare_data(data)

    # Compute Mann-Whitney U statistics dynamically
    u_stroke, p_stroke = mannwhitneyu(stroke_pos, stroke_neg, alternative='two-sided')
    u_bl,     p_bl     = mannwhitneyu(bl_beh_pos, bl_beh_neg, alternative='two-sided')
    u_ac,     p_ac     = mannwhitneyu(ac_beh_pos, ac_beh_neg, alternative='two-sided')

    # Build figure
    fig, axes = plt.subplots(1, 3, figsize=(8.5, 4.5))

    strip_box_panel(
        axes[0],
        data_pos=stroke_pos,
        data_neg=stroke_neg,
        ylabel='Infarct Volume (mm³)',
        panel_label='A',
        p_val=p_stroke,
    )

    strip_box_panel(
        axes[1],
        data_pos=bl_beh_pos,
        data_neg=bl_beh_neg,
        ylabel='Forelimb Asymmetry Score (%)',
        panel_label='B',
        p_val=p_bl,
        zero_line=True,
    )

    strip_box_panel(
        axes[2],
        data_pos=ac_beh_pos,
        data_neg=ac_beh_neg,
        ylabel='Forelimb Asymmetry Score (%)',
        panel_label='C',
        p_val=p_ac,
        zero_line=True,
    )

    plt.tight_layout(pad=0.5)

    save_figure(fig, 'Figure1_InjurySeverity', OUTPUT_DIR)
    plt.close(fig)
    print('Figure 1 complete.')


if __name__ == '__main__':
    main()


