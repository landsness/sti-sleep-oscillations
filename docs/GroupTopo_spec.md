# Spec: Group-average SO-power topo maps (MATLAB)

**Goal.** Add MATLAB code that produces two group-average cortical topo figures of slow-oscillation (SO) power expressed as **% of baseline**, at the **Acute (24 h)** and **Week 1** timepoints, with the lesional/perilesional ROIs and the 1 mm photothrombotic (PT) spot overlaid. These are for the Frontiers revision (Reviewer 1) and are meant to sit beside Figure 4A, so the intensity scale must read as **% Baseline (baseline = 100%)**, matching 4A's y-axis.

Prefer extending the existing script **`SO_GroupTopoplots_AverageFromRawMaps_SavePNGs_v2_arn.m`** (it already averages from the raw maps and saves PNGs). Reuse its map-loading and masking where possible.

---

## Inputs (already exist)

- **Per-mouse maps:** `Mouse_<ID>_SO_RawMaps.mat`, each containing `so_map_bsl`, `so_map_acute`, `so_map_wk1` — all `128×128 double`, **linear** SO-band (0.1–1 Hz) power, already registered to the common bregma-centered atlas frame. Zeros/NaN outside the brain.
- **ROI table:** `ROI_Summary_V2.xlsx`, one row per mouse. Columns used: `Mouse`, `ROI1_CentroidX`, `ROI1_CentroidY`, `ROI2_CentroidX`, `ROI2_CentroidY`, `ROI1_Pixels`, `ROI2_Pixels`. Centroids are in MATLAB image coordinates (x = column, y = row, 1-indexed), from `regionprops`. ROI1 = lesional, ROI2 = perilesional.
- **Mice:** IDs `1–17` and `19–24` (skip 18; exclude the `2DBSI`/`4DBSI` pilots). n = 23.

## Calibration (fixed)

- Registered window = `[-5 5 -5 5]` mm over 128 px ⇒ **mmpp = 10/128 = 0.0781 mm/px**.
- Bregma at image center ≈ `(64.5, 64.5)`; midline column ≈ 64.5.
- 1 mm PT spot ⇒ radius `rL = 0.5/mmpp ≈ 6.4 px`.

## Metric

Per mouse, per pixel, **% of baseline** (NOT % change):

```
pct_acute = 100 * so_map_acute ./ so_map_bsl;
pct_wk1   = 100 * so_map_wk1   ./ so_map_bsl;
```

Baseline = 100% by definition, so only Acute and Week 1 panels are drawn.

Brain mask per mouse: `brain = isfinite(so_map_bsl) & (so_map_bsl > 0.05*max(so_map_bsl(:)));`
Set `pct_*` to `NaN` outside `brain`.

---

## Figure A — whole cortex (atlas-registered)

Straight pixelwise average in the atlas frame; keeps both hemispheres and true anatomy.

1. Stack `pct_acute` and `pct_wk1` across the 23 mice into `128×128×23` arrays.
2. **Coverage mask:** fraction of mice with a valid (non-NaN) pixel; display only where coverage ≥ 0.80.
3. `G_acute = mean(...,3,'omitnan')`, `G_wk1 = mean(...,3,'omitnan')`; blank pixels below coverage.
4. Overlays (single representative set): mean ROI1 centroid = mean over mice of `(ROI1_CentroidX, ROI1_CentroidY)`; mean ROI2 centroid likewise. Circle radii `r1 = sqrt(mean(ROI1_Pixels)/pi)`, `r2 = sqrt(mean(ROI2_Pixels)/pi)`. Draw the PT spot ring (radius `rL`) **centered on the mean ROI1 centroid**. Draw bregma `+` and dashed midline.

## Figure B — lesion-centered (ipsilesional)

Align each mouse on its own ROI1 before averaging, so ROI1 sits at the center of a sharp concentric field.

1. For each mouse, translate its `pct_*` map (and its brain mask) so the ROI1 centroid moves to image center:
   `tx = 64.5 - ROI1_CentroidX; ty = 64.5 - ROI1_CentroidY;`
   `shifted = imtranslate(pct_map, [tx ty], 'FillValues', NaN);`
   (Shift the mask with nearest-neighbor / `'FillValues',0` then threshold.)
2. Average across mice with the same ≥80% coverage rule.
3. Overlays: ROI1 at center `(64.5,64.5)`; ROI2 at center **plus the mean (ROI2−ROI1) offset** across mice; PT spot ring centered at `(64.5,64.5)`. Add a 1 mm scale bar. No bregma/midline (frame is lesion-relative, not anatomical).

---

## Shared styling (match Fig 4A)

- Two side-by-side panels: **Acute (24 h)**, **Week 1**. Panel letter **E**.
- `imagesc(G)`; `axis image off`; set `AlphaData = ~isnan(G)` so out-of-brain is transparent/white. `set(gca,'Color','w')`.
- **Colormap:** blue sequential, dark = suppressed. Use `brewermap(256,'Blues')` reversed (`flipud`) if available, else `flipud(bone)` or a hand-made blue ramp. **`caxis([10 100])`.**
- Shared colorbar labeled **`SO Power (% Baseline)`**, ticks `[20 40 60 80 100]`, with an over-range arrow if any pixel > 100.
- ROI overlays as line circles (parametric `cos/sin`), black with a thin white outline/halo so they read on the blue field: **ROI1 solid, ROI2 dashed, PT spot dotted.** Legend: "ROI1 (lesional)", "ROI2 (perilesional)", "1 mm PT spot".
- Draw the ipsilesional (stroke) hemisphere on the **left** of the image and label **L / R** on Figure A. Optional light brain outline (boundary of the coverage mask).

## Outputs

- Save each figure as PNG (300 dpi) and vector (PDF/EPS) to the results folder, e.g. `GroupTopo_WholeCortex_pctBaseline.(png|pdf)` and `GroupTopo_LesionCentered_pctBaseline.(png|pdf)`.
- Print a sanity readout: value at ROI1 center and ROI2 center for each timepoint. Expected ≈ ROI1 14% (acute) / 41% (wk1); ROI2 27% / 54% — these reproduce the manuscript's 17.2% / 30.8% lesional/perilesional numbers and confirm the frame.

## Notes / gotchas

- MATLAB is 1-indexed; centroids from `regionprops` are already in that convention — do **not** add/subtract 1 when overlaying with `imagesc`.
- Use `'omitnan'` everywhere; guard divide-by-zero (mask requires `so_map_bsl > 0`).
- Keep the ROI masks out of this — masks were never saved; only centroids + pixel counts exist, so ROIs are drawn as equivalent-area circles.
- The reference Python implementation produced both figures from exactly these inputs; the numbers above are the target to match.
