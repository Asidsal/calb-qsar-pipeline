# AUDIT RECORD for `docking_scores_enrichment.csv` - automated reconstruction of the producing cell by the platform lineage tool.
# Absolute local paths were replaced by P(<filename>); nothing else was edited. This file has NOT been re-run.
# Inputs recorded: chemotype_assignments.csv, conformal_predictor.pkl
from pathlib import Path
ROOT = Path(__file__).resolve().parents[2]
DATA_DIR = ROOT / "data"
def P(name):
    """Locate a file by name in data/, results/ or models/ (intermediate files not shipped are noted in INDEX.md)."""
    for d in ("data", "results", "models"):
        p = ROOT / d / name
        if p.exists():
            return str(p)
    return str(ROOT / "results" / name)

import pandas as pd
import numpy as np
from sklearn.metrics import roc_auc_score

ds = pd.read_csv("dock/docking_scores.csv").dropna(subset=['affinity'])
ds['role'] = np.where(ds.tag.str.startswith('active'), 'active', 'decoy')
ds['y'] = (ds['role'] == 'active').astype(int)
ds.to_csv("docking_scores_enrichment.csv", index=False)