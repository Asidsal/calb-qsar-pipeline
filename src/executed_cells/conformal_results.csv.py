# AUDIT RECORD for `conformal_results.csv` - automated reconstruction of the producing cell by the platform lineage tool.
# Absolute local paths were replaced by P(<filename>); nothing else was edited. This file has NOT been re-run.
# Inputs recorded: chemotype_assignments.csv, X_train_ecfp4.npy, X_test_ecfp4.npy, X_train_desc.csv, X_test_desc.csv, y_train.npy, y_test.npy, train.csv, test.csv
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

import numpy as np, pandas as pd, json, warnings
warnings.filterwarnings('ignore')
from sklearn.model_selection import train_test_split
from xgboost import XGBClassifier

SEED = 42

Xtr_fp = np.load(P("X_train_ecfp4.npy"))
Xte_fp = np.load(P("X_test_ecfp4.npy"))
Xtr_d = pd.read_csv(P("X_train_desc.csv"))
Xte_d = pd.read_csv(P("X_test_desc.csv"))
ytr = np.load(P("y_train.npy"))
yte = np.load(P("y_test.npy"))

Xtr = np.hstack([Xtr_fp, Xtr_d.values])
Xte = np.hstack([Xte_fp, Xte_d.values])

def mk():
    sp = (ytr == 0).sum() / ytr.sum()
    return XGBClassifier(n_estimators=500, learning_rate=0.05, max_depth=6, subsample=0.8,
        colsample_bytree=0.8, scale_pos_weight=sp, eval_metric='logloss', random_state=SEED, n_jobs=-1)

# three-way split so an exchangeable test set exists
idx = np.arange(len(ytr))
i_rest, i_rand = train_test_split(idx, test_size=0.20, stratify=ytr, random_state=SEED)
i_prop, i_cal = train_test_split(i_rest, test_size=0.25, stratify=ytr[i_rest], random_state=SEED)
print("proper", len(i_prop), "cal", len(i_cal), "random-test", len(i_rand), "scaffold-test", len(yte))
m = mk(); m.fit(Xtr[i_prop], ytr[i_prop])

def nc(P, lab): return 1.0 - P[:, lab]
def build_cal(groups, a, y):
    return {(g, l): np.sort(a[(groups == g) & (y == l)]) for g in np.unique(groups) for l in (0, 1)}
def pvals(P, groups, cal):
    pv = np.ones((len(P), 2))
    for l in (0, 1):
        a = nc(P, l)
        for g in np.unique(groups):
            mk_ = groups == g; ref = cal.get((g, l))
            if ref is None or len(ref) == 0: continue
            pv[mk_, l] = (len(ref) - np.searchsorted(ref, a[mk_], side='left') + 1) / (len(ref) + 1)
    return pv

P_cal = m.predict_proba(Xtr[i_cal]); y_cal = ytr[i_cal]
a_cal = np.array([1 - P_cal[i, y_cal[i]] for i in range(len(y_cal))])
g0 = np.zeros(len(y_cal), dtype=int)
cal_std = build_cal(g0, a_cal, y_cal)

P_rand = m.predict_proba(Xtr[i_rand]); y_rand = ytr[i_rand]
P_scaf = m.predict_proba(Xte)
pv_rand = pvals(P_rand, np.zeros(len(P_rand), dtype=int), cal_std)
pv_scaf = pvals(P_scaf, np.zeros(len(P_scaf), dtype=int), cal_std)

def summarise(pv, y, eps):
    reg = pv > eps                                  # prediction region
    err = ~reg[np.arange(len(y)), y]               # true label excluded -> error
    single = reg.sum(1) == 1; empty = reg.sum(1) == 0; both = reg.sum(1) == 2
    out = dict(eps=eps, error=err.mean(), singleton=single.mean(), empty=empty.mean(), both=both.mean())
    for l in (0, 1):
        s = y == l; out[f'error_class{l}'] = (~reg[np.arange(len(y)), y])[s].mean()
    if single.sum():
        lab = reg[single].argmax(1); out['acc_singleton'] = (lab == y[single]).mean()
    else: out['acc_singleton'] = np.nan
    return out

eps_grid = np.round(np.arange(0.01, 0.305, 0.01), 3)
rows = []
for e in eps_grid:
    r1 = summarise(pv_rand, y_rand, e); r1['split'] = 'random (exchangeable)'; rows.append(r1)
    r2 = summarise(pv_scaf, yte, e);    r2['split'] = 'scaffold (shifted)';    rows.append(r2)
cp = pd.DataFrame(rows)
cp.to_csv("conformal_results.csv", index=False)
for s in cp.split.unique():
    d = cp[cp.split == s]
    print(f"\n{s}: max |empirical-nominal| = {(d.error - d.eps).abs().max():.4f}")
    print(d[d.eps.isin([0.05, 0.10, 0.20])][['eps', 'error', 'error_class0', 'error_class1', 'singleton', 'empty', 'both', 'acc_singleton']].round(4).to_string(index=False))
np.save("pv_scaf.npy", pv_scaf); np.save("pv_rand.npy", pv_rand); np.save("i_rand.npy", i_rand)