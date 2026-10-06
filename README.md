# Slow Oscillations Track Acute Stroke Injury but Not Functional Recovery: analysis and figure code

Lee J, Jadav AA, Landsness EC. Frontiers in Neurology (in revision).

This repository regenerates every statistic, table and figure in the manuscript from the data in `data/`.

## Quick start

```bash
pip install -r requirements.txt
python make_all.py            # about 4 minutes; everything is written to results/
python make_all.py --no-tiff  # skip the Frontiers TIFF export
```

Every script can also be run on its own, e.g. `python figures/fig4_recovery_dissociation.py`.

## Manuscript item → script

| Manuscript item | Script | Output (under `results/`) |
|---|---|---|
| Figure 1 (A histology; B–D infarct volume and forelimb asymmetry by STI) | `figures/fig1_injury_severity.py` (A: `figures/panels/fig1a_histology.py`) | `figures/*/Figure1_InjurySeverity.*` |
| Figure 2 (A lesion-centred map; B traces; C spectra; D–I acute SO) | `figures/fig2_acute_so.py` (A–C: `figures/panels/fig2{a,b,c}_*.py`) | `figures/*/Figure2_AcuteSO.*` |
| Figure 3 | `figures/fig3_bilateral_covariation.py` | `figures/*/Figure3_BilateralSOCovariation.*` |
| Figure 4 | `figures/fig4_recovery_dissociation.py` | `figures/*/Figure4_RecoveryDissociation.*` |
| Figure 5 | `figures/fig5_baseline_li.py` | `figures/*/Figure5_BaselineLI.*` |
| Supplementary Figure 1 (infarct depth) | `figures/supplfig1_infarct_depth.py` via `analysis/infarct_depth.py` | `figures/*/SupplFig1_InfarctDepth.*`, `tables/infarct_depth_fractions.csv` |
| Supplementary Figure 2 | `figures/supplfig2_infarct_vs_behavior.py` | `figures/*/SupplFig2_InfarctVsBehavior.*` |
| Supplementary Figure 3 | `figures/supplfig3_infarct_vs_acute_so.py` | `figures/*/SupplFig3_InfarctVsAcuteSO.*` |
| Supplementary Figure 4 | `figures/supplfig4_sti_subgroups.py` | `figures/*/SupplFig4_STISubgroups.*` |
| Supplementary Figure 5 | `figures/supplfig5_acute_li.py` | `figures/*/SupplFig5_AcuteLI.*` |
| Supplementary Figure 6 (week-1 SO vs behavior; replaces former Table 1, R1.16) | `figures/supplfig6_wk1_so_vs_behavior.py` | `figures/*/SupplFig6_Wk1SO_vs_Behavior.*` |
| Supplementary Figure 7 | `figures/supplfig7_baseline_so_null.py` | `figures/*/SupplFig7_BaselineSO_Null.*` |
| Supplementary Figure 8 | `figures/supplfig8_baseline_li_peri_null.py` | `figures/*/SupplFig8_BaselineLI_Peri_Null.*` |
| Supplementary Table 1 (n = 23 vs n = 25) | `figures/suppltable1_sensitivity.py` | `tables/SupplTable1_Sensitivity.csv`, `figures/*/SupplTable1_Sensitivity.*` |
| Reviewer analyses R1.11, R1.13, R1.17, R2.1 | `analysis/step1_remaining_analyses.py` | `stats/step1_results.json` |
| Supplementary Table 2 (ρ with 95% CI, p and FDR q for all 34 correlations; R2.2, R2 minor 2) | `analysis/fdr.py` | `tables/SupplTable2_FDR.csv` |
| Every in-text statistic, recomputed and checked against the manuscript, with its effect size and 95% CI (ρ/r: Fisher z; U/W: rank-biserial r, bootstrap) | `analysis/stats_ledger.py` | `stats/stats_ledger.csv` |
| Frontiers submission TIFFs | `figures/export_frontiers.py` | `frontiers/*.tif` |

## Layout

```
config.py         every path, the week-1 exclusion list and the atlas scale, in one place
make_all.py       runs the whole pipeline
data/             all inputs (see data/README.md for provenance and Data@Becker locations)
analysis/         common.py (data loading, partial Spearman, Steiger, effect sizes / CIs), checks, statistics
figures/          one script per manuscript figure; panels/ builds the image panels of Figs 1-2
upstream/         MATLAB that produced the inputs from the raw imaging (not run by make_all.py)
docs/             figure build specifications and captions
archive/          retired scripts, kept for provenance
results/          generated outputs (not version-controlled)
```

## Data provenance (summary)

```
raw WFCI (Data@Becker)
  └─ upstream/matlab/LeeWrapper2Cam_el.m ........... preprocessing -> *-fc1-Power.mat, *-dataFluor.mat
       ├─ SO_PowerAnalysis_Arnav_v2.m .............. data/rawmaps/Mouse_<ID>_SO_RawMaps.mat   (Fig 2A)
       │    └─ SO_ROI_Quantification_mirrored_v2.m . data/ROI_Summary_V2.xlsx (hand-drawn ROIs)
       │         └─ copied by hand ................. data/WFCI_AI_Summary.xlsx (all statistics)
       ├─ Extract_ROI_Spectra.m .................... data/derived/ROI_Spectra_export.mat      (Fig 2C)
       └─ Extract_ROI_Traces.m ..................... data/derived/ROI_Traces_export.mat       (Fig 2B)
histology (Data@Becker)
  ├─ ImageJ infarct areas .......................... data/Stroke_Quantification.xlsx -> infarct volume
  ├─ ImageJ depth tracings (J. Lee) ................ data/depth_traces/*_RoiSet.zip            (Suppl Fig 1)
  └─ representative sections ....................... data/histology/*.tif                     (Fig 1A)
behavior (cylinder test, Data@Becker) .............. data/WFCI_AI_Summary.xlsx (behavior_score)
```

`analysis/check_inputs.py` confirms that the hand-curated `WFCI_AI_Summary.xlsx` matches its sources: all 300 SO-power values match `ROI_Summary_V2.xlsx`, and all 25 infarct volumes match `Stroke_Quantification.xlsx`.

## Data availability

Raw and processed data: Digital Commons Data@Becker, https://doi.org/10.17632/8n5h7v2965.
