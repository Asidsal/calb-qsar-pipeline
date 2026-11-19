# AUDIT RECORD for `conformal_chemotype_validity.csv` - automated reconstruction of the producing cell by the platform lineage tool.
# Absolute local paths were replaced by P(<filename>); nothing else was edited. This file has NOT been re-run.
# Inputs recorded: chemotype_assignments.csv, cal_idx.npy, prop_idx.npy
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
from sklearn.model_selection import train_test_split
from xgboost import XGBClassifier
from rdkit import Chem
from rdkit.Chem import Descriptors, rdMolDescriptors

# Chemotype assignment helpers
PATTERNS = [
    ("Triazole",        ["[CX4]n1cncn1", "[CX4]n1nccn1"]),
    ("Imidazole",       ["[CX4]n1ccnc1"]),
    ("Fluoropyrimidine", ["Nc1nc(=O)[nH]cc1F", "Fc1c[nX2,nX3]c(=O)[nX2,nX3]c1"]),
    ("Allylamine",      ["[NX3;!$(NC=O)](C)CC=C"]),
    ("Thiocarbamate",   ["[NX3]C(=S)O", "[NX3]C(=O)S"]),
    ("Morpholine",      ["C1COCCN1"]),
    ("Benzimidazole",   ["c1ccc2[nH]cnc2c1"]),
]
COMPILED = []
for n, p in PATTERNS:
    ms = [Chem.MolFromSmarts(s) for s in p]
    bad = [s for s, m in zip(p, ms) if m is None]
    assert not bad, (n, bad)
    COMPILED.append((n, ms))

def chemotype(smi):
    m = Chem.MolFromSmiles(smi)
    if m is None: return "Unparsed"
    mw = Descriptors.MolWt(m)
    ri = m.GetRingInfo()
    rings = [len(r) for r in ri.AtomRings()]
    if mw > 800 and rdMolDescriptors.CalcNumAmideBonds(m) >= 5: return "Echinocandin"
    if rings and max(rings) >= 20 and mw > 600: return "Polyene"
    for n, ms in COMPILED:
        if any(m.HasSubstructMatch(q) for q in ms): return n
    return "Other"

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

tr = pd.read_csv(R + r"\train.csv")
te = pd.read_csv(R + r"\test.csv")

chem = pd.read_csv(P("chemotype_assignments.csv"))[['std_smiles', 'chemotype']].drop_duplicates('std_smiles')
tr = tr.merge(chem, on='std_smiles', how='left')
te = te.merge(chem, on='std_smiles', how='left')

def mk():
    sp = (ytr == 0).sum() / ytr.sum()
    return XGBClassifier(n_estimators=500, learning_rate=0.05, max_depth=6, subsample=0.8,
        colsample_bytree=0.8, scale_pos_weight=sp, eval_metric='logloss', random_state=SEED, n_jobs=-1)

# three-way split
idx = np.arange(len(ytr))
i_rest, i_rand = train_test_split(idx, test_size=0.20, stratify=ytr, random_state=SEED)
i_prop, i_cal = train_test_split(i_rest, test_size=0.25, stratify=ytr[i_rest], random_state=SEED)

m = mk()
m.fit(Xtr[i_prop], ytr[i_prop])

def nc(P, lab): return 1.0 - P[:, lab]

def build_cal(groups, a, y):
    return {(g, l): np.sort(a[(groups == g) & (y == l)]) for g in np.unique(groups) for l in (0, 1)}

def pvals(P, groups, cal):
    pv = np.ones((len(P), 2))
    for l in (0, 1):
        a = nc(P, l)
        for g in np.unique(groups):
            mk_ = groups == g
            ref = cal.get((g, l))
            if ref is None or len(ref) == 0: continue
            pv[mk_, l] = (len(ref) - np.searchsorted(ref, a[mk_], side='left') + 1) / (len(ref) + 1)
    return pv

P_cal = m.predict_proba(Xtr[i_cal])
y_cal = ytr[i_cal]
a_cal = np.array([1 - P_cal[i, y_cal[i]] for i in range(len(y_cal))])
g0 = np.zeros(len(y_cal), dtype=int)
cal_std = build_cal(g0, a_cal, y_cal)

P_scaf = m.predict_proba(Xte)
pv_scaf = pvals(P_scaf, np.zeros(len(P_scaf), dtype=int), cal_std)

# chemotype-conditional Mondrian
MIN_GRP = 40
chem_tr = tr.chemotype.values
chem_te = te.chemotype.values
g_cal_raw = chem_tr[i_cal]
counts = pd.Series(g_cal_raw).value_counts()
keep = set(counts[counts >= MIN_GRP].index)

def grp(a): return np.array([x if x in keep else 'POOLED' for x in a])

g_cal = grp(g_cal_raw)
g_te = grp(chem_te)
cal_chem = build_cal(g_cal, a_cal, y_cal)
for l in (0, 1):
    if len(cal_chem.get(('POOLED', l), [])) == 0:
        cal_chem[('POOLED', l)] = np.sort(a_cal[y_cal == l])
pv_chem = pvals(P_scaf, g_te, cal_chem)

def per_chemotype_error(pv, y, chem, eps):
    reg = pv > eps
    err = ~reg[np.arange(len(y)), y]
    out = []
    for c in pd.unique(chem):
        s = chem == c
        if s.sum() < 20: continue
        out.append(dict(chemotype=c, n=int(s.sum()), error=float(err[s].mean())))
    return pd.DataFrame(out)

EPS = 0.10
std_e = per_chemotype_error(pv_scaf, yte, chem_te, EPS).rename(columns={'error': 'error_class_conditional'})
chm_e = per_chemotype_error(pv_chem, yte, chem_te, EPS).rename(columns={'error': 'error_chemotype_conditional'})
cmp_ = std_e.merge(chm_e[['chemotype', 'error_chemotype_conditional']], on='chemotype').sort_values('n', ascending=False)
cmp_['nominal'] = EPS
cmp_.to_csv("conformal_chemotype_validity.csv", index=False)