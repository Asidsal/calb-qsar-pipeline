# AUDIT RECORD for `candidates_with_docking.csv` - automated reconstruction of the producing cell by the platform lineage tool.
# Absolute local paths were replaced by P(<filename>); nothing else was edited. This file has NOT been re-run.
# Inputs recorded: -
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

import os, glob, subprocess, re, time, csv
import numpy as np, pandas as pd
import warnings
warnings.filterwarnings('ignore')
from rdkit import Chem, DataStructs, RDLogger
from rdkit.Chem import BRICS, Descriptors, QED, rdMolDescriptors, AllChem, Crippen, RDConfig
from rdkit.Chem.Scaffolds import MurckoScaffold
from rdkit.Chem.FilterCatalog import FilterCatalog, FilterCatalogParams
from rdkit.ML.Descriptors import MoleculeDescriptors
from rdkit.DataStructs import ConvertToNumpyArray
import sys, pickle, random
RDLogger.DisableLog('rdApp.*')
sys.path.append(os.path.join(RDConfig.RDContribDir, 'SA_Score'))
import sascorer

R = str(DATA_DIR)

# Load candidates with docking scores
cand = pd.read_csv("docking_candidates.csv")
cs = pd.read_csv(os.path.join("dock", "candidate_scores.csv")).dropna(subset=['affinity'])
m = cand.merge(cs, on='cand_id', how='left')
m.to_csv("candidates_with_docking.csv", index=False)