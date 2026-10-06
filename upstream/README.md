# upstream/ — MATLAB that produced the inputs in data/

`make_all.py` does not run these scripts. They need the full preprocessed imaging set (hundreds of GB, on Data@Becker) and, for the ROI step, interactive drawing. They are included so the whole chain from raw imaging to `data/` is documented. The paths at the top of each script point to the lab's original locations (`Z:\...`, `Box\DBSI Directory\...`); edit them to your local copy of the deposit.

## matlab/ (the pipeline, in order)

| Step | Script | Input | Output |
|---|---|---|---|
| 1 | `LeeWrapper2Cam_el.m` + `LeeWrapper2Cam/` helpers (`procFluor`, `gsr`, `PowerAnalysis`, `qcRaw`, …), run lists via `parseRunsJM.m` | raw WFCI | `PreProcessed_Data/<ID>/<session>/*-fc1-Power.mat`, `*-dataFluor.mat`, `*-LandmarksAndMask.mat` |
| 2 | `SO_PowerAnalysis_Arnav_v2.m` | `*-fc1-Power.mat` (`whole_spectra_map`, contrast 4 = corrected calcium; mean over 0.1–1 Hz) | `Mouse_<ID>_SO_RawMaps.mat` → `data/rawmaps/` |
| 3 | `SO_ROI_Quantification_mirrored_v2.m` + `DrawROI_arn_v2.m`, `autoCaxis.m` (interactive) | RawMaps | one row per mouse in `ROI_Summary_V2.xlsx` → `data/` |
| 4 | `Extract_ROI_Spectra.m` | Power.mat + ROI_Summary_V2 | `ROI_Spectra_export.mat` → `data/derived/` (Fig 2C) |
| 4 | `Extract_ROI_Traces.m` | dataFluor.mat + ROI_Summary_V2 | `ROI_Traces_export.mat` → `data/derived/` (Fig 2B) |
| — | `SO_GroupTopoplots_AverageFromRawMaps_SavePNGs_v2_arn.m` | RawMaps | exploratory group topoplots. The manuscript's Figure 2A is produced in Python by `figures/panels/fig2a_lesion_topo.py`. |

## supporting/ (reviewer-response checks; no figure depends on them)

| Script | What it checks |
|---|---|
| `SO_Contrast_SO_Maps.m`, `SO_Contrast_SO_Maps_direct.m`, `SO_Contrast_SO_Maps_direct_1.m`, `so_contrast_maps.py` | SO maps across imaging contrasts (corrected vs uncorrected calcium, hemoglobin) |
| `SO_Uncorrected_vs_Corrected.m` | Effect of hemodynamic correction on SO power |
| `SO_Artifact_BandCheck.m` | Band-power check of the bright-core artifact (R1 major critique) |
