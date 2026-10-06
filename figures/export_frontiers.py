"""Export all manuscript figures as Frontiers in Neurology-compliant TIFF files.

Frontiers specs:
  - TIFF, LZW compression, RGB 8-bit, 300 dpi
  - Two-column max: 7.09 in (180 mm); single-column: 3.35 in (85 mm)
  - Minimum font size: 8 pt; fonts/markers scaled to match JNeurosci proportions

Output directory: results/frontiers/
    Figure1.tif ... Figure5.tif, SupplementaryFigure1.tif ... SupplementaryFigure8.tif

This is the ONLY place the script-stem -> submission-filename mapping lives;
if figures are renumbered, edit STEM_MAP and SCRIPTS below.
Image panels in Figures 1-2 (results/panels/) are built first by make_all.py.
"""

import sys
import io
import importlib
from pathlib import Path

FIGURES_DIR = Path(__file__).resolve().parent
sys.path[:0] = [str(FIGURES_DIR), str(FIGURES_DIR.parent), str(FIGURES_DIR.parent / 'analysis')]
import config  # noqa: E402

import matplotlib.pyplot as plt
from PIL import Image
import figure_style

OUT_DIR = config.FRONTIERS
OUT_DIR.mkdir(parents=True, exist_ok=True)

# ---------------------------------------------------------------------------
# Output stem -> Frontiers filename (same mapping as JNeurosci)
# ---------------------------------------------------------------------------
STEM_MAP = {
    'Figure1_InjurySeverity':            'Figure1',
    'Figure2_AcuteSO':                   'Figure2',
    'Figure3_BilateralSOCovariation':    'Figure3',
    'Figure4_RecoveryDissociation':      'Figure4',
    'Figure5_BaselineLI':                'Figure5',
    'SupplFig1_InfarctDepth':            'SupplementaryFigure1',
    'SupplFig2_InfarctVsBehavior':       'SupplementaryFigure2',
    'SupplFig3_InfarctVsAcuteSO':        'SupplementaryFigure3',
    'SupplFig4_STISubgroups':            'SupplementaryFigure4',
    'SupplFig5_AcuteLI':                 'SupplementaryFigure5',
    'SupplFig6_Wk1SO_vs_Behavior':       'SupplementaryFigure6',
    'SupplFig7_BaselineSO_Null':         'SupplementaryFigure7',
    'SupplFig8_BaselineLI_Peri_Null':    'SupplementaryFigure8',
}

# ---------------------------------------------------------------------------
# Target figsize per script
# Original sizes: Fig1=(8.5,4.5) Fig2=(15,9) Fig3=(10,9) Fig4=(18,4.8)
#                 Fig5=(10,4.8)  S1-S2,S4=(10,4.5) S3=(10,13) S5=(5,4.8)
# Two-column: 7.09 in (180 mm); single-column: 3.35 in (85 mm)
# Figure 3: reduced nominal to 6.6 in to stay within 180 mm after tight bbox
# ---------------------------------------------------------------------------
MAX_H = 9.5
TWO_COL = 7.09
ONE_COL = 3.35

def _h(orig_w, orig_h, tgt_w):
    return min(tgt_w * orig_h / orig_w, MAX_H)

TARGET_FIGSIZE = {
    'fig1_injury_severity':              (TWO_COL, 9.1),                       # histology row + 3 panels
    'fig2_acute_so':                     (6.95, 5.8),                          # image row + 2 x 3 panels
    'fig3_bilateral_covariation':        (6.6,     _h(10.0, 9.0, 6.6)),
    'fig4_recovery_dissociation':        (TWO_COL, 2.8),                       # proportional is too flat
    'fig5_baseline_li':                  (TWO_COL, _h(10.0, 4.8, TWO_COL)),
    'supplfig1_infarct_depth':           (TWO_COL, _h(10.0, 3.4, TWO_COL)),
    'supplfig2_infarct_vs_behavior':     (TWO_COL, _h(10.0, 4.5, TWO_COL)),
    'supplfig3_infarct_vs_acute_so':     (TWO_COL, _h(9.2,  8.4, TWO_COL)),
    'supplfig4_sti_subgroups':           (TWO_COL, _h(10.0, 4.5, TWO_COL)),
    'supplfig5_acute_li':                (TWO_COL, _h(10.0, 13.0, TWO_COL)),   # capped at MAX_H
    'supplfig6_wk1_so_vs_behavior':      (TWO_COL, _h(9.2,  8.4, TWO_COL)),
    'supplfig7_baseline_so_null':        (TWO_COL, _h(9.2,  8.4, TWO_COL)),   # 2 x 2 (rev. 2026-10-05)
    'supplfig8_baseline_li_peri_null':   (ONE_COL, _h(5.0,  4.8, ONE_COL)),    # single column
}

DPI = 300
results = []

# ---------------------------------------------------------------------------
# Style constants scaled for 7" print width — same as JNeurosci export.
# Frontiers minimum is 8 pt; bump tick/stats/legend to 8 (from JNeurosci's 6).
# ---------------------------------------------------------------------------
FR_FONT_SIZES = {
    'panel_label': 8,
    'title':       8,
    'axis_label':  8,
    'tick':        8,
    'stats':       8,
    'legend':      8,
    'sig_star':    8,
}
FR_LINE_WIDTHS = {
    'regression':  0.75,
    'identity':    0.5,
    'axis':        0.5,
    'bracket':     0.6,
    'errorbar':    0.8,
    'median':      1.2,
    'box_whisker': 0.6,
}
FR_MARKER_SIZES = {
    'scatter':  20,
    'errorbar':  4,
    'strip':    15,
}

orig_font_sizes   = dict(figure_style.FONT_SIZES)
orig_line_widths  = dict(figure_style.LINE_WIDTHS)
orig_marker_sizes = dict(figure_style.MARKER_SIZES)

figure_style.FONT_SIZES.update(FR_FONT_SIZES)
figure_style.LINE_WIDTHS.update(FR_LINE_WIDTHS)
figure_style.MARKER_SIZES.update(FR_MARKER_SIZES)

_orig_subplots   = plt.subplots
_current_figsize = [None]

def _patched_subplots(*args, **kwargs):
    if _current_figsize[0] is not None:
        kwargs['figsize'] = _current_figsize[0]
    return _orig_subplots(*args, **kwargs)

plt.subplots = _patched_subplots

# ---------------------------------------------------------------------------
# Patch save_figure to also write Frontiers TIFF
# ---------------------------------------------------------------------------
_orig_save_figure = figure_style.save_figure

def _patched_save_figure(fig, output_stem, outdir, formats=('svg', 'pdf', 'png')):
    # Write ONLY the Frontiers TIFF; results/figures/ keeps the standard-size
    # outputs written by make_all.py (the old version overwrote them).
    _export_tiff(fig, output_stem)
    return {}

def _export_tiff(fig, stem):
    jn_name = STEM_MAP.get(stem)
    if not jn_name:
        print(f'  [WARN] No Frontiers mapping for "{stem}" -- skipped')
        return

    buf = io.BytesIO()
    fig.savefig(buf, format='png', dpi=DPI, bbox_inches='tight',
                pad_inches=0.05, facecolor='white')
    buf.seek(0)
    img = Image.open(buf).convert('RGB')
    max_px = int(180 / 25.4 * DPI)          # 2125 px = 180 mm at 300 dpi (Frontiers two-column max)
    if img.size[0] > max_px:                # tight bbox can overshoot by a few px; scale down, keep dpi
        img = img.resize((max_px, round(img.size[1] * max_px / img.size[0])), Image.LANCZOS)

    tif_path = OUT_DIR / f'{jn_name}.tif'
    img.save(str(tif_path), format='TIFF', compression='tiff_lzw', dpi=(DPI, DPI))

    w_px, h_px = img.size
    w_mm = w_px / DPI * 25.4
    size_mb = tif_path.stat().st_size / 1e6

    results.append({
        'file':  f'{jn_name}.tif',
        'source': stem,
        'px':    f'{w_px}x{h_px}',
        'mm':    f'{w_mm:.0f}',
        'dpi':   DPI,
        'mb':    f'{size_mb:.1f}',
    })
    ok = 'OK' if w_mm <= 180.05 else 'WARN: exceeds 180 mm'
    print(f'  -> {jn_name}.tif  [{w_px}x{h_px} px, {w_mm:.0f} mm wide]  {ok}')

figure_style.save_figure = _patched_save_figure

# ---------------------------------------------------------------------------
# Run each script
# ---------------------------------------------------------------------------
SCRIPTS = list(TARGET_FIGSIZE)

for script_name in SCRIPTS:
    print(f'\n=== {script_name}  target: {TARGET_FIGSIZE[script_name]} ===')
    _current_figsize[0] = TARGET_FIGSIZE[script_name]
    mod = importlib.import_module(script_name)
    mod.main()
    plt.close('all')

plt.subplots = _orig_subplots
figure_style.FONT_SIZES.update(orig_font_sizes)
figure_style.LINE_WIDTHS.update(orig_line_widths)
figure_style.MARKER_SIZES.update(orig_marker_sizes)

# ---------------------------------------------------------------------------
# Summary
# ---------------------------------------------------------------------------
print('\n' + '=' * 72)
print(f'{"File":<16} {"Source stem":<38} {"Pixels":<14} {"mm wide":<8} {"MB"}')
print('-' * 72)
for r in results:
    print(f'{r["file"]:<16} {r["source"]:<38} {r["px"]:<14} {r["mm"]:<8} {r["mb"]}')
print('=' * 72)
print(f'\nAll files written to: {OUT_DIR}')
