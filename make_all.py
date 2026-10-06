"""make_all.py -- regenerate every statistic, table and figure in the manuscript.

    python make_all.py            # full run (about 4 minutes)
    python make_all.py --no-tiff  # skip the Frontiers TIFF export

Steps
  1. analysis/check_inputs.py     curated table vs machine-generated sources
  2. analysis/infarct_depth.py    infarct depth fractions (Suppl Fig 1, R1 major)
     analysis/fdr.py              BH-FDR over the correlation family -> Suppl Table 2 (figures read q from here)
  3. figures/panels/*.py          image panels for Figures 1A, 2A, 2B, 2C
  4. figures/fig*.py, supplfig*.py, suppltable1_sensitivity.py
  5. analysis/step1_remaining_analyses.py  reviewer analyses R1.11, R1.13, R1.17, R2.1
     analysis/stats_ledger.py     every test reported in the manuscript, recomputed, with effect size + 95% CI
  6. figures/export_frontiers.py  submission TIFFs (results/frontiers/)
All outputs go to results/ (see README.md).
"""
import os
import runpy
import sys
import time
from pathlib import Path

ROOT = Path(__file__).resolve().parent
os.environ.setdefault('MPLBACKEND', 'Agg')

STEPS = [
    'analysis/check_inputs.py',
    'analysis/infarct_depth.py',
    'analysis/fdr.py',
    'figures/panels/fig1a_histology.py',
    'figures/panels/fig2a_lesion_topo.py',
    'figures/panels/fig2b_lesional_traces.py',
    'figures/panels/fig2c_roi_spectra.py',
    'figures/fig1_injury_severity.py',
    'figures/fig2_acute_so.py',
    'figures/fig3_bilateral_covariation.py',
    'figures/fig4_recovery_dissociation.py',
    'figures/fig5_baseline_li.py',
    'figures/supplfig1_infarct_depth.py',
    'figures/supplfig2_infarct_vs_behavior.py',
    'figures/supplfig3_infarct_vs_acute_so.py',
    'figures/supplfig4_sti_subgroups.py',
    'figures/supplfig5_acute_li.py',
    'figures/supplfig6_wk1_so_vs_behavior.py',
    'figures/supplfig7_baseline_so_null.py',
    'figures/supplfig8_baseline_li_peri_null.py',
    'figures/suppltable1_sensitivity.py',
    'analysis/step1_remaining_analyses.py',
    'analysis/stats_ledger.py',
]
if '--no-tiff' not in sys.argv:
    STEPS.append('figures/export_frontiers.py')


def run(rel):
    path = ROOT / rel
    old_argv, old_path = sys.argv[:], sys.path[:]
    sys.argv = [str(path)]
    sys.path.insert(0, str(path.parent))
    try:
        runpy.run_path(str(path), run_name='__main__')
    finally:
        sys.argv, sys.path[:] = old_argv, old_path
        import matplotlib.pyplot as plt
        plt.close('all')


if __name__ == '__main__':
    t0 = time.time()
    for rel in STEPS:
        print(f'\n{"=" * 72}\n>>> {rel}\n{"=" * 72}', flush=True)
        run(rel)
    print(f'\nAll steps completed in {time.time() - t0:.0f} s. Outputs: {ROOT / "results"}')
