# Figure Generation Pipeline

## Quick Start — Regenerate All Figures

```bash
cd figures/

# Main figures
python make_figure_1.py
python make_figure_2.py
python make_figure_3.py
python make_figure_4.py
python make_figure_5.py

# Extended Data figures
python make_suppl_s1.py
python make_suppl_s2.py
python make_suppl_s3.py
python make_suppl_s4.py
python make_suppl_s5.py
```

## JNeurosci Submission Export

To regenerate all ten submission-ready TIFF files at once:

```bash
cd figures/
python export_jneurosci.py
```

Outputs → `figures_jneurosci/` (created automatically):

| File | Source script | Content |
|---|---|---|
| `Figure1.tif` | `make_figure_1.py` | Injury severity |
| `Figure2.tif` | `make_figure_2.py` | Acute SO suppression |
| `Figure3.tif` | `make_figure_3.py` | Bilateral SO covariation |
| `Figure4.tif` | `make_figure_4.py` | SO recovery and dissociation |
| `Figure5.tif` | `make_figure_5.py` | Baseline LI predicts recovery |
| `Figure1-1.tif` | `make_suppl_s1.py` | Extended Data: infarct vs behavior |
| `Figure2-1.tif` | `make_suppl_s2.py` | Extended Data: STI subgroup robustness |
| `Figure2-2.tif` | `make_suppl_s3.py` | Extended Data: acute LI vs outcomes |
| `Figure5-1.tif` | `make_suppl_s4.py` | Extended Data: baseline SO null |
| `Figure5-2.tif` | `make_suppl_s5.py` | Extended Data: baseline LI ROI2 null |

**TIFF specs:** RGB 8-bit, LZW compression, 600 dpi, 7.0" wide (3.4" for Figure5-2),
fonts 6–7 pt, markers and line weights scaled for print size.

The export script re-runs each generating script at the JNeurosci print dimensions
(fonts, markers, and line weights scaled; figsize overridden per figure). It does not
convert existing PNGs. Each run also refreshes `figures/final/` as a side effect.

Figure 4 target height is set to 2.8" (overriding proportional scaling, which produces
a figure too flat for a 1×4 layout at 7" width).

---

## Output

Running individual scripts saves to `figures/final/`:
- `svg/`  — vector, editable (`svg.fonttype=none` — text as text, not paths)
- `pdf/`  — for manuscript submission
- `png/`  — 600 dpi preview

---

## Requirements

Python 3.x with: `matplotlib`, `pandas`, `scipy`, `numpy`, `openpyxl`, `Pillow`

```bash
pip install matplotlib pandas scipy numpy openpyxl pillow
```

---

## Data File

`WFCI AI Summary.xlsx` (sheet: `"AI Summary "` — note trailing space) must be in the project
root (`BiomniJune2/`). Scripts resolve the path automatically via `Path(__file__).parent.parent`.

- n = 25 animals (STI+ n=14, STI− n=11)
- Outliers **DBSI_M02** and **DBSI_M04** excluded from **Week 1 SO power analyses only** (n=23)
- All behavioral correlations use n=25

---

## Shared Style Library

`figure_style.py` — call `set_publication_style()` at the top of any script. All constants
(`COLORS`, `FONT_SIZES`, `LINE_WIDTHS`, `MARKER_SIZES`, `EXPORT_SETTINGS`) and helper
functions (`scatter_with_regression`, `strip_box_panel`, `add_panel_label`, `format_stats`,
`load_data`, `save_figure`, …) are defined here.

---

## Figure Map

| Script | Output stem | Panels | Content |
|---|---|---|---|
| `make_figure_1.py` | `Figure1_InjurySeverity` | A–C (1×3) | Injury severity: infarct volume, baseline behavior, acute behavior by STI group (strip+box) |
| `make_figure_2.py` | `Figure2_AcuteSO` | A–F (2×3) | **Row 0:** A=lateralization/gradient strip+box (all animals, no legend — colors match Figure 3 row/col labels); B=STI-group strip+box; C=ROI1 ipsi vs acute behavior. **Row 1:** D=ROI2 ipsi vs acute behavior; E=ROI1 contra vs acute behavior (null, stats lower right); F=ROI2 contra vs acute behavior (null) |
| `make_figure_3.py` | `Figure3_BilateralSOCovariation` | A–D (2×2 matrix) | **Baseline top row, acute bottom row; ROI1 left col, ROI2 right col.** A=BL ROI1; B=BL ROI2; C=Acute ROI1; D=Acute ROI2. Row/column labels added. |
| `make_figure_4.py` | `Figure4_RecoveryDissociation` | A–D (1×4) | A=SO power recovery trajectory (4-line % baseline, n=23); B=bilateral SO covariation trajectory (Pearson r, n=23); C=STI fractional SO recovery strip+box (null); D=STI week-1 behavior strip+box (significant) |
| `make_figure_5.py` | `Figure5_BaselineLI` | A–B (1×2) | A=Baseline LI ROI1 vs Week-1 behavior (Spearman + partial ρ); B=Baseline LI ROI1 vs infarct volume (independence check, null) |
| `make_suppl_s1.py` | `SupplFig_S1_InfarctVsBehavior` | A–B (1×2) | Infarct size vs acute behavior (A) and week-1 behavior (B) |
| `make_suppl_s2.py` | `SupplFig_S2_STISubgroups` | A–B (1×2) | STI-subgroup robustness: acute ROI1 ipsi SO vs acute behavior within STI+ (A) and STI− (B) |
| `make_suppl_s3.py` | `SupplFig_S3_AcuteLI` | A–F (3×2) | Acute LI across all outcome measures. **Row 0:** A=LI ROI1 vs acute beh; B=LI ROI2 vs acute beh. **Row 1:** C=LI ROI1 vs week-1 beh; D=LI ROI2 vs week-1 beh. **Row 2:** E=LI ROI1 vs infarct; F=LI ROI2 vs infarct |
| `make_suppl_s4.py` | `SupplFig_S4_BaselineSO_Null` | A–B (1×2) | Baseline absolute SO does not predict week-1 recovery: ROI1 (A), ROI2 (B) — null result |
| `make_suppl_s5.py` | `SupplFig_S5_BaselineLI_ROI2_Null` | A (1×1) | Baseline LI ROI2 vs week-1 behavior — null result (scatter) |

---

## Manuscript Citation Notes

### Figure 3 panel order (important)
Figure 3 uses a **matrix layout**: baseline is the top row, acute is the bottom row; ROI1 is the
left column, ROI2 is the right column.

| | ROI1 (left) | ROI2 (right) |
|---|---|---|
| **Baseline (top)** | A | B |
| **Acute 24h (bottom)** | C | D |

In-text references:
- Both baseline values → "Figure 3A,B"
- Acute ROI1 → "Figure 3C"
- Acute ROI2 → "Figure 3D"

---

## Archive

Retired scripts and pre-reorganization outputs are preserved in `figures/archive/`.
