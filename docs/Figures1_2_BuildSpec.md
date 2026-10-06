# Figures 1 & 2 — Build Spec (revision handoff)

**Status:** current revised design for Figures 1 and 2 of the STI/SO manuscript, as of the R1-major-critique revision. NOT final publication quality — statistics (p-values, ρ, partial correlations, depth numbers) may change while other reviewer critiques (e.g., multiple-comparison correction) are addressed. This spec exists so the figures can be regenerated cleanly if numbers or inputs change.

This document is written for Claude Code (or any operator) to rebuild the figures from the scripts and data below.

---

## 0. Global conventions

- **Output formats:** each figure/panel rendered to PNG (300 dpi) + vector PDF + SVG.
- **Fonts:** sans-serif; panel letters bold. `svg.fonttype='none'`, `pdf.fonttype=42`.
- **Palette (Okabe–Ito, colorblind-safe):**
  - Timepoints: Baseline = black `#000000`; Acute = vermillion `#D55E00`; Week 1 = blue `#0072B2`.
  - Uncorrected trace = orange `#E69F00`; Corrected trace = black `#000000`.
  - STI groups (existing boxplot/scatter panels): STI+ = red, STI− = green (as in original manuscript figures).
  - SO band shading = `#F5D400`, alpha 0.13; band = 0.1–1 Hz.
- **Spectra axes:** log frequency, log power; SO band shaded; ~8 Hz cardiac peak may be annotated.
- **Spatial scale (topoplots/ROIs):** registered atlas window [-5 5 -5 5] mm over 128 px → 78.1 µm/px, bregma-centered.
- **ROI geometry:** ROI1 (lesional) ~0.58 mm dia; ROI2 (perilesional) ~0.37 mm dia, centroid ~0.80 mm from ROI1 (medial + slightly anterior); PT laser spot = 1.0 mm dia (0.5 mm radius), centered on ROI1.

## 1. Data sources (on user's machine / Box)

- **Preprocessed WFCI:** `...\WF_Calcium_RepositoryREDONE\PreProcessed_Data\<mouse>\{bsl,acute,oneweek}\`
  - `*-fc1-Power.mat` → `whole_spectra_map[128,128,freq,contrast]`, `hz`; contrast 4 = hemodynamically corrected calcium.
  - `*-fc1-dataFluor.mat` → `xform_datafluor` (uncorrected), `xform_datafluorCorr` (corrected), `runInfo.samplingRate` (≈20 Hz).
  - `*-LandmarksAndMask.mat` → `xform_isbrain`.
- **ROI centroids/sizes:** `...\SO_Analysis\ROI_Summary_V2.xlsx` (`ROI1_CentroidX/Y`, `ROI1_Pixels`, `ROI2_*`). Contralateral homolog = ROI1 mirrored across midline (x → 129−x).
- **Infarct areas per slice:** `...\histology\Stroke Quanitification.xlsx` (per-mouse `Ms#` sheets: `Slice`, `Stroke Infarct Area (micron^2)`; exclude the `Total` and `Additional` rows).
- **CV histology TIFFs:** `...\histology\cresyl_violet\M##\WFCI_Ms#_<slice>.tif` (folders zero-padded, filenames not).
- **STI status + stroke volume:** `WFCI AI Summary.xlsx` (per-mouse `secondary_thalamic_injury` 0/1, `stroke_size` mm³).
- **Cohort:** SO analyses use n = 23 (mice 1–17, 19–24; no 18). Mouse 1 = "no stroke" but included in the n = 23. Behavior–infarct correlations use n = 25.

## 2. Scripts (in this folder)

| Script | Lang | Input | Output |
|---|---|---|---|
| `Extract_ROI_Spectra.m` | MATLAB | Power.mat tree + ROI_Summary_V2 | `ROI_Spectra_export.mat` (ROI-avg PSD, all mice/sessions) |
| `Extract_ROI_Traces.m` | MATLAB | dataFluor.mat tree + ROI_Summary_V2 | `ROI_Traces_export.mat` (ROI-avg uncorr+corr traces) |
| `plot_roi_spectra.py` | Python | `ROI_Spectra_export.mat` | Fig 2C spectra (3 ROIs, bsl/acute[/wk1]) |
| `plot_roi_traces.py` | Python | `ROI_Traces_export.mat` | Fig 2B traces (Ms8 lesional, bsl/acute[/wk1]) |
| `make_figure1A_montage.py` | Python | 4 section images | Fig 1A histology 2×2 montage |
| `GroupTopo_MATLAB_spec_for_Codex.md` | spec | group lesion-centered SO maps | Fig 2A topoplot source |

MATLAB extractions run where the big data live (user's machine); they emit small export files. Python plotting runs on the exports. Set the `INPUT`/path constants at the top of each script.

---

## 3. FIGURE 1 — Injury severity (supports §3.1)

Panels:
- **A — Histology montage (NEW).** 2×2: columns STI+/STI−, rows cresyl violet (cortical infarct) / GFAP (ipsilesional thalamus). Build: `make_figure1A_montage.py`.
  - **Representative animals (LOCKED):**
    - Cresyl violet STI+ = **M11, peak slice 8** → raw `histology\cresyl_violet\M11\WFCI_Ms11_8.tif`.
    - Cresyl violet STI− = **M22, peak slice 13** → raw `histology\cresyl_violet\M22\WFCI_Ms22_13.tif`.
    - GFAP STI+ = **M15** → `histology\IHC\M15\WFCI_Ms15_2.tif`.
    - GFAP STI− = **M02** → `histology\IHC\M02\Ms2_1.tif`.
  - **Source images:** use the **RAW** TIFFs above (no burned-in text/scale bars). Earlier annotated `CV_STIpos_*/CV_STIneg_*.png` derivatives are deprecated — they carried a burned-in mouse/volume label (top-left) and an "1 mm (approx)" bar (bottom-left) that surfaced in the top-right/bottom-right after the flip. The montage draws its own scale bars, so raw inputs are required.
  - **Orientation (LOCKED):** all four panels are **horizontally mirrored** (`PIL.ImageOps.mirror`) because the strokes are **left-hemisphere**. After the flip the affected (left) hemisphere sits on the **viewer's left** in every panel, and the GFAP thalamic gliosis is ipsilateral to (directly below) the cortical infarct. Secondary thalamic injury is always ipsilateral to the cortical lesion (thalamocortical projections are ipsilateral), so CV infarct and GFAP gliosis must be on the same side.
  - **Arrows (STI+ GFAP only):** two white `-|>` arrows with a black outline stroke, pointing up-left onto the two right-thalamic gliotic foci of the *mirrored* image. Coordinates are in **percent of image W×H** on the flipped Ms15_2 (1920×1440): arrow 1 head (26, 44) / tail (21, 49.5); arrow 2 head (16, 58) / tail (11.5, 63.5). Style: `arrowstyle='-|>'`, `color='white'`, `lw=2.2`, `mutation_scale=24`, `shrinkA=0`, `shrinkB=1`, `path_effects=[withStroke(linewidth=4, foreground='black')]`. (These are the pre-flip coordinates head1 74/44, tail1 79/49.5, head2 84/58, tail2 88.5/63.5 mirrored via x→100−x; arrows were shortened to half length and thinned so arrow 1's tail no longer overlaps arrow 2's target.)
  - **Scale bars — TODO for final:** the montage currently draws a fixed decorative "1 mm" bar per tile (150 px on the 900 px tile); CV and GFAP are at different magnifications, so calibrate each bar to the true µm/px before the final figure.
- **B — Infarct volume by STI** (was panel A). Source: `WFCI AI Summary` stroke_size × STI. Stat (current): Mann–Whitney U = 151, p < 0.001.
- **C — Baseline forelimb asymmetry by STI** (was B). Current: U = 82, p = 0.805.
- **D — Acute forelimb asymmetry by STI** (was C). Current: U = 106, p = 0.119.

Panels B–D come from the **original manuscript figure-generation code** (boxplots + points, STI+ red / STI− green). This spec adds histology as A and shifts B/C/D. Supplementary Figure 1 (infarct vs behavior: ρ = 0.636 acute, ρ = 0.502 wk1) unchanged.

**May change:** none of the Fig 1 stats are multiple-comparison-sensitive in the current plan, but re-verify B–D if group definitions change.

## 4. FIGURE 2 — Acute SO suppression (supports §3.2)

Order (top row new, then existing):
- **A — Acute/baseline SO topoplot (NEW).** Group, lesion-centered, % baseline, with ROI1/ROI2/laser overlay. Build from group lesion-centered SO maps per `GroupTopo_MATLAB_spec_for_Codex.md`; overlay rings to scale (78.1 µm/px). **TODO:** regenerate as a clean single ACUTE panel (own colorbar, 1 mm scale bar) rather than the current crop `panel_acute_topo.png`.
- **B — Representative lesional traces (NEW).** Ms8, baseline vs acute, uncorrected (orange) vs corrected (black) ΔF/F, 60 s, y shared. Build: `plot_roi_traces.py` (`MOUSE=8`, sessions bsl+acute). **Ms8 locked** (highest SO SNR among near-median mice 8/9/14/24).
- **C — ROI power spectra (NEW).** Lesional/perilesional/contralateral, baseline vs acute, group mean ± SEM (n = 23), log–log, SO band shaded. Build: `plot_roi_spectra.py` with `SESS=['bsl','acute']`. (A bsl+acute+wk1 variant exists for other uses.)
- **D — Acute SO power by region×hemisphere** (% baseline) (was A). Current stats: all Wilcoxon p < 0.001; ipsi<contra (lesion W=3, peri W=8, both p<0.001); peri>lesion (30.8% vs 17.2%, W=1, p<0.001).
- **E — Same, STI+ vs STI−** (was B). Lesion U=38 p=0.035; peri U=45 p=0.085; contra lesion U=78 p=0.978, peri U=93 p=0.396.
- **F–I — Acute SO vs forelimb asymmetry** (were C–F): F ipsi-lesion (ρ=−0.718, p<0.001), G ipsi-peri (ρ=−0.600, p=0.002), H contra-lesion (ρ=−0.047, p=0.824), I contra-peri (ρ=−0.079, p=0.707); n=25. Partial-ρ (controlling infarct size): lesion −0.448 (p=0.025), peri −0.248 (p=0.232).

Panels D–I come from the **original manuscript figure code**; this spec adds A–C and relabels the originals D–I.

**LIKELY TO CHANGE:** the correlation and partial-correlation statistics in **E/F–I** (and Supplementary Fig 2) are the ones most exposed to multiple-comparison correction and other-reviewer edits. When they change, re-run the original correlation code and update panel annotations; the new A–C panels (traces/spectra/topoplot) are descriptive and will not change unless ROI definitions or the cohort change.

---

## 5. Regeneration workflow

1. If underlying data or ROI definitions change: re-run `Extract_ROI_Spectra.m` and `Extract_ROI_Traces.m` (MATLAB, on the data machine) → refresh the two export `.mat` files.
2. `python plot_roi_spectra.py` → Fig 2C; `python plot_roi_traces.py` → Fig 2B.
3. Regenerate Fig 2A topoplot from the group maps (MATLAB spec) with ROI/laser overlay.
4. `python make_figure1A_montage.py` → Fig 1A (points at the raw CV + GFAP TIFFs; mirrors all panels; draws the STI+ arrows).
5. Re-run original manuscript code for Fig 1B–D and Fig 2D–I if their stats changed.
6. Assemble panels (Illustrator/InDesign or a layout script) into final Figure 1 and Figure 2; apply panel letters A… in the orders above.

## 6. Locked decisions
- Trace representative: **Ms8**. Histology representatives (CV): **STI+ = M11 (slice 8), STI− = M22 (slice 13)**; GFAP: **STI+ = M15 (WFCI_Ms15_2), STI− = M02 (Ms2_1)**.
- Histology montage uses **RAW** TIFFs (no annotated derivatives) and is **horizontally mirrored** (left-hemisphere strokes → affected hemisphere on viewer's left; GFAP gliosis ipsilateral/below the CV infarct).
- Figure 2 panel order: topoplot (A) → traces (B) → spectra (C) → boxplots (D) → STI (E) → behavior (F–I).
- Figure 1 panel order: histology (A) → infarct volume (B) → baseline behavior (C) → acute behavior (D).

## 7. Open placeholders
- Fig 1A: DONE — final raw CV + GFAP images in place, mirrored to left-affected orientation, arrows on STI+ thalamic gliosis. Remaining: calibrate scale bars (see §3 panel A TODO).
- Depth analysis: DONE — formal coronal lesion-incidence heatmap from Jake's ImageJ traces (`SupplFig_DepthHeatmap.{png,pdf,svg}`, STI+ N=14 / STI− N=11). Reported on the N=25 cohort to match the manuscript: 24 of 25 had a delineable infarct, 19 of these spared deep cortex (depth fraction ≤ 0.8; median 0.57). (Do NOT cite a bare N=22 anywhere — it confuses the reviewer against the stated 25-mouse cohort.) **Locked as Supplementary Figure 1.**
- Fig 2A topoplot: replace crop with clean regenerated acute panel.
- Supplementary figure numbering LOCKED: depth = S1 (cited first in §3.1); prior figures shift up one — old S1→S2 (infarct vs. behavior), S2→S3 (STI subgroups), S3→S4 (acute LI), S4→S5 (baseline SO null), S5→S6 (baseline LI ROI2 null). Manuscript body citations updated (tracked); the separate Supplementary Material file and `SupplFig_S#_*` image files still need renumbering to match.
