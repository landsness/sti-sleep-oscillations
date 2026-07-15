"""Export all ten manuscript figures as Frontiers of Neurology-compliant TIFF files.

Frontiers specs:
  - TIFF, LZW compression, RGB 8-bit, 300 dpi
  - Two-column max: 7.09 in (180 mm); single-column: 3.35 in (85 mm)
  - Minimum font size: 8 pt; fonts/markers scaled to match JNeurosci proportions

Output directory: figures_frontiers/
    Figure1.tif  Figure2.tif  Figure3.tif  Figure4.tif  Figure5.tif
    Figure1-1.tif  Figure2-1.tif  Figure2-2.tif  Figure5-1.tif  Figure5-2.tif
"""

import sys
import io
import importlib
from pathlib import Path

FIGURES_DIR = Path(__file__).parent
sys.path.insert(0, str(FIGURES_DIR))

import matplotlib.pyplot as plt
from PIL import Image
import figure_style

OUT_DIR = FIGURES_DIR.parent / 'figures_frontiers'
OUT_DIR.mkdir(exist_ok=True)

# ---------------------------------------------------------------------------
# Output stem -> Frontiers filename (same mapping as JNeurosci)
# ---------------------------------------------------------------------------
STEM_MAP = {
    'Figure1_InjurySeverity':           'Figure1',
    'Figure2_AcuteSO':                  'Figure2',
    'Figure3_BilateralSOCovariation':   'Figure3',
    'Figure4_RecoveryDissociation':     'Figure4',
    'Figure5_BaselineLI':               'Figure5',
    'SupplFig_S1_InfarctVsBehavior':    'Figure1-1',
    'SupplFig_S2_STISubgroups':         'Figure2-1',
    'SupplFig_S3_AcuteLI':             'Figure2-2',
    'SupplFig_S4_BaselineSO_Null':      'Figure5-1',
    'SupplFig_S5_BaselineLI_ROI2_Null': 'Figure5-2',
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
    'make_figure_1': (TWO_COL, _h(8.5,  4.5, TWO_COL)),   # 7.09 x 3.75
    'make_figure_2': (TWO_COL, 5.0),                        # extra height for 6 crowded panels
    'make_figure_3': (6.6,     _h(10.0, 9.0, 6.6)),        # 6.6 nominal -> ~180 mm after tight bbox
    'make_figure_4': (TWO_COL, 2.8),                        # override: proportional too flat
    'make_figure_5': (TWO_COL, _h(10.0, 4.8, TWO_COL)),   # 7.09 x 3.40
    'make_suppl_s1': (TWO_COL, _h(10.0, 4.5, TWO_COL)),   # 7.09 x 3.19
    'make_suppl_s2': (TWO_COL, _h(10.0, 4.5, TWO_COL)),   # 7.09 x 3.19
    'make_suppl_s3': (TWO_COL, _h(10.0,13.0, TWO_COL)),   # 7.09 x 9.22 (capped at 9.5)
    'make_suppl_s4': (TWO_COL, _h(10.0, 4.5, TWO_COL)),   # 7.09 x 3.19
    'make_suppl_s5': (ONE_COL, _h(5.0,  4.8, ONE_COL)),   # 3.35 x 3.22 (single column)
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
    result = _orig_save_figure(fig, output_stem, outdir, formats)
    _export_tiff(fig, output_stem)
    return result

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
    ok = 'OK' if w_mm <= 181 else 'WARN: exceeds 180 mm'
    print(f'  -> {jn_name}.tif  [{w_px}x{h_px} px, {w_mm:.0f} mm wide]  {ok}')

figure_style.save_figure = _patched_save_figure

# ---------------------------------------------------------------------------
# Run each script
# ---------------------------------------------------------------------------
SCRIPTS = [
    'make_figure_1',
    'make_figure_2',
    'make_figure_3',
    'make_figure_4',
    'make_figure_5',
    'make_suppl_s1',
    'make_suppl_s2',
    'make_suppl_s3',
    'make_suppl_s4',
    'make_suppl_s5',
]

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
