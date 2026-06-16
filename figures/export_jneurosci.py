"""Export all ten manuscript figures as JNeurosci-compliant TIFF files.

Each figure is re-rendered at its JNeurosci print size with scaled fonts,
markers, and line weights — not converted from existing PNGs.

Output directory: BiomniJune2/figures_jneurosci/
    Figure1.tif  Figure2.tif  Figure3.tif  Figure4.tif  Figure5.tif
    Figure1-1.tif  Figure2-1.tif  Figure2-2.tif  Figure5-1.tif  Figure5-2.tif

JNeurosci specs:
  - TIFF, LZW compression, RGB 8-bit, 600 dpi
  - Double-column max: 7.0 in wide; single-column: 3.4 in (Figure5-2 only)
  - Height scaled proportionally from original figsize, capped at 9.5 in
  - Fonts 6-8 pt, line weights 0.5-0.8 pt, markers scaled to print size
"""

import sys
import io
import importlib
from pathlib import Path

FIGURES_DIR = Path(__file__).parent
sys.path.insert(0, str(FIGURES_DIR))

import matplotlib.pyplot as plt
from PIL import Image
import figure_style  # import before patching

OUT_DIR = FIGURES_DIR.parent / 'figures_jneurosci'
OUT_DIR.mkdir(exist_ok=True)

# ---------------------------------------------------------------------------
# Output stem -> JNeurosci filename
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
# Target figsize per script (width scaled to JNeurosci, height proportional)
# Original sizes: Fig1=(8.5,4.5) Fig2=(15,9) Fig3=(10,9) Fig4=(18,4.8)
#                 Fig5=(10,4.8)  S1-S2,S4=(10,4.5) S3=(10,13) S5=(5,4.8)
# ---------------------------------------------------------------------------
MAX_H = 9.5

def _h(orig_w, orig_h, tgt_w):
    return min(tgt_w * orig_h / orig_w, MAX_H)

TARGET_FIGSIZE = {
    'make_figure_1': (7.0, _h(8.5,  4.5, 7.0)),   # 7.0 x 3.71
    'make_figure_2': (7.0, _h(15.0, 9.0, 7.0)),   # 7.0 x 4.20
    'make_figure_3': (6.5, _h(10.0, 9.0, 6.5)),   # 6.5 nominal -> ~17.56 cm actual (tight bbox adds ~0.4"); confirmed compliant
    'make_figure_4': (7.0, 2.8),                   # override: proportional scaling too flat
    'make_figure_5': (7.0, _h(10.0, 4.8, 7.0)),   # 7.0 x 3.36
    'make_suppl_s1': (7.0, _h(10.0, 4.5, 7.0)),   # 7.0 x 3.15
    'make_suppl_s2': (7.0, _h(10.0, 4.5, 7.0)),   # 7.0 x 3.15
    'make_suppl_s3': (7.0, _h(10.0,13.0, 7.0)),   # 7.0 x 9.10 (capped at 9.5)
    'make_suppl_s4': (7.0, _h(10.0, 4.5, 7.0)),   # 7.0 x 3.15
    'make_suppl_s5': (3.4, _h(5.0,  4.8, 3.4)),   # 3.4 x 3.26
}

# ---------------------------------------------------------------------------
# JNeurosci-appropriate style constants (at final print size)
# ---------------------------------------------------------------------------
JN_FONT_SIZES = {
    'panel_label': 8,    # bold A B C
    'title':       7,
    'axis_label':  7,
    'tick':        6,
    'stats':       6,
    'legend':      6,
    'sig_star':    7,
}

JN_LINE_WIDTHS = {
    'regression':  0.75,
    'identity':    0.5,
    'axis':        0.5,
    'bracket':     0.6,
    'errorbar':    0.8,
    'median':      1.2,
    'box_whisker': 0.6,
}

JN_MARKER_SIZES = {
    'scatter':  20,
    'errorbar':  4,
    'strip':    15,
}

DPI = 600
results = []

# ---------------------------------------------------------------------------
# Patch figure_style dicts IN-PLACE so all already-imported scripts see changes
# ---------------------------------------------------------------------------
orig_font_sizes   = dict(figure_style.FONT_SIZES)
orig_line_widths  = dict(figure_style.LINE_WIDTHS)
orig_marker_sizes = dict(figure_style.MARKER_SIZES)

figure_style.FONT_SIZES.update(JN_FONT_SIZES)
figure_style.LINE_WIDTHS.update(JN_LINE_WIDTHS)
figure_style.MARKER_SIZES.update(JN_MARKER_SIZES)

# ---------------------------------------------------------------------------
# Patch plt.subplots to override figsize
# ---------------------------------------------------------------------------
_orig_subplots   = plt.subplots
_current_figsize = [None]  # mutable cell for closure

def _patched_subplots(*args, **kwargs):
    if _current_figsize[0] is not None:
        kwargs['figsize'] = _current_figsize[0]
    return _orig_subplots(*args, **kwargs)

plt.subplots = _patched_subplots

# ---------------------------------------------------------------------------
# Patch save_figure to also write JNeurosci TIFF
# ---------------------------------------------------------------------------
_orig_save_figure = figure_style.save_figure

def _patched_save_figure(fig, output_stem, outdir, formats=('svg', 'pdf', 'png')):
    result = _orig_save_figure(fig, output_stem, outdir, formats)
    _export_tiff(fig, output_stem)
    return result

def _export_tiff(fig, stem):
    jn_name = STEM_MAP.get(stem)
    if not jn_name:
        print(f'  [WARN] No JNeurosci mapping for "{stem}" -- skipped')
        return

    # Render to lossless PNG buffer, convert to RGB, save as LZW TIFF
    buf = io.BytesIO()
    fig.savefig(buf, format='png', dpi=DPI, bbox_inches='tight',
                pad_inches=0.05, facecolor='white')
    buf.seek(0)
    img = Image.open(buf).convert('RGB')

    tif_path = OUT_DIR / f'{jn_name}.tif'
    img.save(str(tif_path), format='TIFF', compression='tiff_lzw', dpi=(DPI, DPI))

    w_px, h_px = img.size
    size_mb = tif_path.stat().st_size / 1e6
    if size_mb > 20:
        print(f'  [WARN] {jn_name}.tif is {size_mb:.1f} MB -- exceeds 20 MB limit')

    results.append({
        'file':  f'{jn_name}.tif',
        'source': stem,
        'px':    f'{w_px}x{h_px}',
        'dpi':   DPI,
        'mode':  'RGB',
        'mb':    f'{size_mb:.1f}',
    })
    print(f'  -> {jn_name}.tif  [{w_px}x{h_px} px, {size_mb:.1f} MB]')

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

# Restore
plt.subplots = _orig_subplots
figure_style.FONT_SIZES.update(orig_font_sizes)
figure_style.LINE_WIDTHS.update(orig_line_widths)
figure_style.MARKER_SIZES.update(orig_marker_sizes)

# ---------------------------------------------------------------------------
# Summary
# ---------------------------------------------------------------------------
print('\n' + '=' * 78)
print(f'{"File":<16} {"Source stem":<38} {"Pixels":<14} {"DPI":<5} {"Mode":<5} {"MB"}')
print('-' * 78)
for r in results:
    print(f'{r["file"]:<16} {r["source"]:<38} {r["px"]:<14} {r["dpi"]:<5} {r["mode"]:<5} {r["mb"]}')
print('=' * 78)
print(f'\nAll files written to: {OUT_DIR}')
