# STI and Sleep Oscillations — Figure Reproduction Code

Code to reproduce all figures from:

> **Acute Slow Oscillation Power as a Biomarker of Injury Severity After Photothrombotic Stroke: Dissociation from Week 1 Functional Recovery**


---

## Repository Contents

```
WFCI AI Summary.xlsx    — source data (n=25 animals)
figures/                — figure generation scripts
  figure_style.py       — shared style library
  make_figure_1.py … make_figure_5.py
  make_suppl_s1.py … make_suppl_s5.py
  export_jneurosci.py   — batch export to JNeurosci-compliant TIFFs
  run_so_recovery_stats.py  — SO recovery statistics (Wilcoxon)
  README.md             — detailed figure map and script documentation
requirements.txt
```

---

## Quickstart

```bash
# 1. Install dependencies
pip install -r requirements.txt

# 2. Generate a single figure (saves SVG/PDF/PNG to figures/final/)
cd figures/
python make_figure_1.py

# 3. Regenerate all figures
python make_figure_1.py && python make_figure_2.py && python make_figure_3.py && \
python make_figure_4.py && python make_figure_5.py && \
python make_suppl_s1.py && python make_suppl_s2.py && python make_suppl_s3.py && \
python make_suppl_s4.py && python make_suppl_s5.py

# 4. Export submission-ready TIFFs (JNeurosci specs: 600 dpi, LZW, RGB)
python export_jneurosci.py
```

Output TIFFs are written to `figures_jneurosci/` (created automatically).

---

## Data

`WFCI AI Summary.xlsx` — sheet `"AI Summary "` (note trailing space in sheet name).

- **n = 25** animals (STI+ n=14, STI− n=11)
- Animals **DBSI_M02** and **DBSI_M04** are excluded from Week 1 SO power analyses only (n=23); all behavioral correlations use n=25
- Column definitions are described in the Methods section of the manuscript

---

## Figure Map

| Script | Manuscript figure | Content |
|---|---|---|
| `make_figure_1.py` | Figure 1 | Injury severity: infarct volume, baseline and acute behavior by STI group |
| `make_figure_2.py` | Figure 2 | Acute SO suppression: lateralization, STI-group differences, behavioral correlations |
| `make_figure_3.py` | Figure 3 | Bilateral SO covariation at baseline and acute timepoints (2×2 matrix) |
| `make_figure_4.py` | Figure 4 | SO recovery trajectory and STI-group dissociation |
| `make_figure_5.py` | Figure 5 | Baseline laterality index predicts week-1 recovery |
| `make_suppl_s1.py` | Extended Data Figure 1-1 | Infarct size vs acute and week-1 behavior |
| `make_suppl_s2.py` | Extended Data Figure 2-1 | STI-subgroup robustness |
| `make_suppl_s3.py` | Extended Data Figure 2-2 | Acute laterality index vs all outcomes |
| `make_suppl_s4.py` | Extended Data Figure 5-1 | Baseline absolute SO: null result |
| `make_suppl_s5.py` | Extended Data Figure 5-2 | Baseline LI ROI2: null result |

See [`figures/README.md`](figures/README.md) for panel-level descriptions and citation notes.

---

## Requirements

Python 3.9+ with packages listed in `requirements.txt`. No GPU or special hardware required.

---

## License

Code: [MIT License](LICENSE)

Data (`WFCI AI Summary.xlsx`): available for research use; please cite the manuscript if used.

---

## Contact

Eric Landsness — landsness@wustl.edu
