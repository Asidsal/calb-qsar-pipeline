# AUDIT RECORD for `cluster_split_repeats.csv` - automated reconstruction of the producing cell by the platform lineage tool.
# Absolute local paths were replaced by P(<filename>); nothing else was edited. This file has NOT been re-run.
# Inputs recorded: pool_meta.csv, y_all.npy
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

import numpy as np
import pandas as pd
import warnings
warnings.filterwarnings('ignore')
from sklearn.metrics import matthews_corrcoef, roc_auc_score, balanced_accuracy_score, confusion_matrix
from xgboost import XGBClassifier

meta = pd.read_csv(P("pool_meta.csv"))
yall = np.load(P("y_all.npy"))

R = str(DATA_DIR)
SEED = 42
Xtr_fp = np.load(R + r"\X_train_ecfp4.npy")
Xte_fp = np.load(R + r"\X_test_ecfp4.npy")
Xtr_d = pd.read_csv(R + r"\X_train_desc.csv")
Xte_d = pd.read_csv(R + r"\X_test_desc.csv")
ytr = np.load(R + r"\y_train.npy")
yte = np.load(R + r"\y_test.npy")
Xtr = np.hstack([Xtr_fp, Xtr_d.values])
Xte = np.hstack([Xte_fp, Xte_d.values])
Xall = np.vstack([Xtr, Xte])

def mk():
    sp = (ytr == 0).sum() / ytr.sum()
    return XGBClassifier(n_estimators=500, learning_rate=0.05, max_depth=6, subsample=0.8,
                         colsample_bytree=0.8, scale_pos_weight=sp, eval_metric='logloss',
                         random_state=SEED, n_jobs=-1)

cl = meta.groupby('cluster').size()
rep = []
for s in range(5):
    r = np.random.default_rng(100 + s)
    order = r.permutation(cl.index.values)
    cum = 0
    tc = set()
    for c in order:
        if cum >= 0.20 * len(meta):
            break
        tc.add(c)
        cum += cl[c]
    tm = meta.cluster.isin(tc).values
    m = mk()
    m.fit(Xall[~tm], yall[~tm])
    yp = m.predict(Xall[tm])
    pr = m.predict_proba(Xall[tm])[:, 1]
    y = yall[tm]
    tn, fp, fn, tp = confusion_matrix(y, yp, labels=[0, 1]).ravel()
    rep.append(dict(seed=100 + s, n_test=int(tm.sum()), active_rate=float(y.mean()),
                    MCC=matthews_corrcoef(y, yp), BA=balanced_accuracy_score(y, yp),
                    AUC=roc_auc_score(y, pr),
                    Sens=tp / (tp + fn), Spec=tn / (tn + fp)))

d = pd.DataFrame(rep)
d.to_csv("cluster_split_repeats.csv", index=False)