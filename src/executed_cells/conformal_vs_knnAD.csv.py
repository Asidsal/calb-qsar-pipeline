# AUDIT RECORD for `conformal_vs_knnAD.csv` - automated reconstruction of the producing cell by the platform lineage tool.
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

import numpy as np, pandas as pd, warnings
warnings.filterwarnings('ignore')
from sklearn.model_selection import train_test_split
from xgboost import XGBClassifier
from rdkit import Chem, DataStructs
from rdkit.Chem import AllChem

R = str(DATA_DIR)
SEED=42

Xtr_fp=np.load(R+r"\X_train_ecfp4.npy"); Xte_fp=np.load(R+r"\X_test_ecfp4.npy")
Xtr_d=pd.read_csv(R+r"\X_train_desc.csv"); Xte_d=pd.read_csv(R+r"\X_test_desc.csv")
ytr=np.load(R+r"\y_train.npy"); yte=np.load(R+r"\y_test.npy")
Xtr=np.hstack([Xtr_fp,Xtr_d.values]); Xte=np.hstack([Xte_fp,Xte_d.values])
tr=pd.read_csv(R+r"\train.csv"); te=pd.read_csv(R+r"\test.csv")

chem=pd.read_csv(P("chemotype_assignments.csv"))[['std_smiles','chemotype']].drop_duplicates('std_smiles')
tr=tr.merge(chem,on='std_smiles',how='left'); te=te.merge(chem,on='std_smiles',how='left')

def mk():
    sp=(ytr==0).sum()/ytr.sum()
    return XGBClassifier(n_estimators=500,learning_rate=0.05,max_depth=6,subsample=0.8,
        colsample_bytree=0.8,scale_pos_weight=sp,eval_metric='logloss',random_state=SEED,n_jobs=-1)

idx=np.arange(len(ytr))
i_rest,i_rand=train_test_split(idx,test_size=0.20,stratify=ytr,random_state=SEED)
i_prop,i_cal=train_test_split(i_rest,test_size=0.25,stratify=ytr[i_rest],random_state=SEED)

m=mk(); m.fit(Xtr[i_prop],ytr[i_prop])

def nc(P,lab): return 1.0-P[:,lab]
def build_cal(groups,a,y):
    return {(g,l):np.sort(a[(groups==g)&(y==l)]) for g in np.unique(groups) for l in (0,1)}
def pvals(P,groups,cal):
    pv=np.ones((len(P),2))
    for l in (0,1):
        a=nc(P,l)
        for g in np.unique(groups):
            mk_=groups==g; ref=cal.get((g,l))
            if ref is None or len(ref)==0: continue
            pv[mk_,l]=(len(ref)-np.searchsorted(ref,a[mk_],side='left')+1)/(len(ref)+1)
    return pv

P_cal=m.predict_proba(Xtr[i_cal]); y_cal=ytr[i_cal]
a_cal=np.array([1-P_cal[i,y_cal[i]] for i in range(len(y_cal))])
g0=np.zeros(len(y_cal),dtype=int)
cal_std=build_cal(g0,a_cal,y_cal)

P_scaf=m.predict_proba(Xte)
pv_scaf=pvals(P_scaf,np.zeros(len(P_scaf),dtype=int),cal_std)

gen=AllChem.GetMorganGenerator(radius=2,fpSize=2048)
def fps(sm):
    o=[]
    for s in sm:
        mol=Chem.MolFromSmiles(s); o.append(gen.GetFingerprint(mol) if mol else None)
    return o
prop_fps=[f for f in fps(tr.std_smiles.values[i_prop]) if f is not None]
te_fps=fps(te.std_smiles.values)
knn5=np.array([np.mean(sorted(DataStructs.BulkTanimotoSimilarity(f,prop_fps),reverse=True)[:5]) if f else 0.0 for f in te_fps])
pred=(P_scaf[:,1]>=0.5).astype(int)
rows=[]
for t in np.arange(0.20,0.86,0.01):
    s=knn5>=t
    if s.sum()<30: continue
    rows.append(dict(method='kNN applicability domain',param=round(float(t),3),
                     retention=float(s.mean()),accuracy=float((pred[s]==yte[s]).mean())))
for e in np.round(np.arange(0.01,0.45,0.01),3):
    reg=pv_scaf>e; single=reg.sum(1)==1
    if single.sum()<30: continue
    lab=reg[single].argmax(1)
    rows.append(dict(method='Mondrian conformal',param=float(e),
                     retention=float(single.mean()),accuracy=float((lab==yte[single]).mean())))
cv=pd.DataFrame(rows); cv.to_csv("conformal_vs_knnAD.csv",index=False)