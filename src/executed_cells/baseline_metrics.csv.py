# AUDIT RECORD for `baseline_metrics.csv` - automated reconstruction of the producing cell by the platform lineage tool.
# Absolute local paths were replaced by P(<filename>); nothing else was edited. This file has NOT been re-run.
# Inputs recorded: repro_svm.json, repro_test_metrics.json, model_results.json
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

import json, numpy as np, pandas as pd

pub=json.load(open(P("model_results.json")))
rep=json.load(open(P("repro_test_metrics.json"))); rep['SVM_ECFP4']=json.load(open(P("repro_svm.json")))
order=['RF_ECFP4','RF_Desc','RF_Combined','SVM_ECFP4','LR_ECFP4','XGB_ECFP4','XGB_Combined']
rows=[]
for m in order:
    for k in ['MCC','AUC','BA','Sens','Spec']:
        rows.append(dict(model=m,metric=k,reported=pub[m]['test_'+k],reproduced=round(float(rep[m][k]),4),
                         delta=round(float(rep[m][k])-pub[m]['test_'+k],4)))
bm=pd.DataFrame(rows)
bm['abs_delta']=bm.delta.abs()
bm.to_csv("baseline_metrics.csv",index=False)