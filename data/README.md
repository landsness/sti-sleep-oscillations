# data/ — inputs and their provenance

Everything the pipeline reads is in this folder. "Deposit" paths refer to the Digital Commons Data@Becker dataset (https://doi.org/10.17632/8n5h7v2965), which mirrors the lab's `DBSI Directory`.

| File | What it is | Produced by | Deposit location |
|---|---|---|---|
| `WFCI_AI_Summary.xlsx` (sheet `AI Summary `, trailing space) | Curated analysis table: one row per mouse × timepoint. Holds infarct volume (`stroke_size`, mm³), STI status, forelimb asymmetry (`behavior_score`) and SO power per ROI/hemisphere. **Every statistic is computed from this table.** | Assembled by hand from the three sources below. Checked by `analysis/check_inputs.py`. | `WFCI AI Summary.xlsx` (top level) |
| `ROI_Summary_V2.xlsx` | Per-mouse ROI1 (lesional) and ROI2 (perilesional) SO power (bsl/acute/wk1 × ipsi/contra), centroids and pixel counts. Contralateral values use the ROI mirrored across the midline (`fliplr`). | `upstream/matlab/SO_ROI_Quantification_mirrored_v2.m` run on the RawMaps. ROIs are hand-drawn polygons (`DrawROI_arn_v2.m`); the masks were not saved, only centroids and pixel counts. | `WFCI/WF_Calcium_RepositoryREDONE/SO_Analysis/ROI_Summary_V2.xlsx` |
| `Stroke_Quantification.xlsx` | Per-slice infarct area (µm²) on cresyl violet sections. Volume = Σ area × 40 µm (the "Total" row). | ImageJ tracing | `histology/Stroke Quanitification.xlsx` |
| `rawmaps/Mouse_<ID>_SO_RawMaps.mat` (25 files) | `so_map_bsl`, `so_map_acute`, `so_map_wk1`: 128 × 128 linear SO-band (0.1–1 Hz) power maps in the bregma-centred atlas frame (10 mm / 128 px = 78.1 µm/px), hemodynamically corrected calcium | `upstream/matlab/SO_PowerAnalysis_Arnav_v2.m` from each session's `*-fc1-Power.mat` | `WFCI/WF_Calcium_RepositoryREDONE/PreProcessed_Data/<ID>/` |
| `derived/ROI_Spectra_export.mat` | ROI-averaged power spectra, all mice and sessions (Fig 2C) | `upstream/matlab/Extract_ROI_Spectra.m` (reads the Power.mat files and ROI_Summary_V2) | derived; regenerate from `PreProcessed_Data` |
| `derived/ROI_Traces_export.mat` | ROI-averaged uncorrected and corrected ΔF/F traces (Fig 2B) | `upstream/matlab/Extract_ROI_Traces.m` (reads the dataFluor.mat files) | derived; regenerate from `PreProcessed_Data` |
| `depth_traces/<mouse>_RoiSet.zip` (24 files) | ImageJ tracings by J. Lee on the peak coronal cresyl violet section, 4 ROIs in order: infarct (polygon), pia, gray/white boundary ("wm"), midline (polylines). Ms1 had no delineable infarct. The Ms23 ROIs are unnamed but in the same order. | ImageJ (manual) | `ImageJ Mice/Ms<N>/` (**not yet deposited**) |
| `histology/*.tif` | Raw sections used in Figure 1A: CV STI+ `WFCI_Ms11_8.tif`, CV STI− `WFCI_Ms22_13.tif`, GFAP STI+ `WFCI_Ms15_2.tif`, GFAP STI− `Ms2_1.tif` | Keyence BZ-X800 | `histology/cresyl_violet/M11/`, `…/M22/`, `histology/IHC/M15/`, `…/M02/` |

## Mouse IDs

| Context | Format |
|---|---|
| Analysis table | `M01`–`M24` (no M18), plus the two pilot animals `DBSI_M02` and `DBSI_M04` |
| RawMaps / ROI summary | `1`–`24`, `2DBSI`, `4DBSI` |
| Depth tracings | `Ms<N>`, `DBSI_Ms2`, `DBSI_Ms4` |

The pilot animals `DBSI_M02` and `DBSI_M04` are excluded from all week-1 SO analyses (`config.WEEK1_SO_EXCLUDED`; see Methods).
