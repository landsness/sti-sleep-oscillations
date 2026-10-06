"""config.py -- every path used by the pipeline, in one place.

All inputs live in data/ (see data/README.md for provenance); all outputs are
written under results/ (git-ignored; regenerate with `python make_all.py`).
"""
from pathlib import Path

ROOT = Path(__file__).resolve().parent

# ---------------------------------------------------------------- inputs
DATA = ROOT / 'data'
XLSX = DATA / 'WFCI_AI_Summary.xlsx'          # curated analysis table (sheet 'AI Summary ')
XLSX_SHEET = 'AI Summary '                    # note trailing space
ROI_SUMMARY = DATA / 'ROI_Summary_V2.xlsx'    # MATLAB ROI output (centroids, pixels, SO power)
RAWMAPS = DATA / 'rawmaps'                    # Mouse_<ID>_SO_RawMaps.mat (128x128 SO-power maps)
DEPTH_TRACES = DATA / 'depth_traces'          # ImageJ RoiSet.zip per mouse (infarct, pia, wm, midline)
HISTOLOGY = DATA / 'histology'                # 4 raw TIFFs used in Figure 1A
DERIVED = DATA / 'derived'                    # small MATLAB exports (ROI spectra / traces)
ROI_SPECTRA = DERIVED / 'ROI_Spectra_export.mat'
ROI_TRACES = DERIVED / 'ROI_Traces_export.mat'

# ---------------------------------------------------------------- outputs
RESULTS = ROOT / 'results'
FIG_OUT = RESULTS / 'figures'                 # {png,pdf,svg}/<stem>.*
PANELS = RESULTS / 'panels'                   # raster/vector panels assembled into Fig 1 and Fig 2
TABLES = RESULTS / 'tables'                   # CSV tables (Suppl Table 1, depth fractions, ...)
STATS = RESULTS / 'stats'                     # statistics ledger / checks
FRONTIERS = RESULTS / 'frontiers'             # submission TIFFs (export_frontiers.py)

# ---------------------------------------------------------------- cohort
# Two pilot animals received a structural MRI under deep isoflurane between the
# acute and week-1 sessions; their week-1 SO is excluded (see Methods).
WEEK1_SO_EXCLUDED = ['DBSI_M02', 'DBSI_M04']
# The same two animals as they are labelled in ROI_Summary_V2.xlsx and in the
# MATLAB exports (data/derived/*.mat): folders Mouse_2DBSI / Mouse_4DBSI.
WEEK1_SO_EXCLUDED_ROI_IDS = ['2DBSI', '4DBSI']

# WFCI atlas frame: [-5 5 -5 5] mm over 128 px, bregma-centred
MM_PER_PX = 10 / 128                          # 0.0781 mm/px
PT_SPOT_RADIUS_MM = 0.5                       # 1-mm photothrombotic spot


def ensure_dirs():
    for d in (FIG_OUT, PANELS, TABLES, STATS, FRONTIERS):
        d.mkdir(parents=True, exist_ok=True)
