# AUDIT RECORD for `generalization_benchmark.csv` - automated reconstruction of the producing cell by the platform lineage tool.
# Absolute local paths were replaced by P(<filename>); nothing else was edited. This file has NOT been re-run.
# Inputs recorded: chemotype_assignments.csv, generalization_sens.csv, _gb_partial.csv
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
from sklearn.metrics import matthews_corrcoef, roc_auc_score, balanced_accuracy_score, confusion_matrix
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

Xall=np.vstack([Xtr,Xte]); yall=np.concatenate([ytr,yte])
meta=pd.concat([tr.assign(orig='train'),te.assign(orig='test')],ignore_index=True)

gen=AllChem.GetMorganGenerator(radius=2,fpSize=2048)
mols=[Chem.MolFromSmiles(s) for s in meta.std_smiles]
fps=[gen.GetFingerprint(m) if m else None for m in mols]

def leader_cluster(fps, cutoff=0.4, seed=42):
    rng=np.random.default_rng(seed)
    order=rng.permutation(len(fps))
    centers=[]; assign=np.full(len(fps),-1)
    for i in order:
        f=fps[i]
        if f is None: continue
        if centers:
            sims=np.array(DataStructs.BulkTanimotoSimilarity(f,[fps[c] for c in centers]))
            j=int(sims.argmax())
            if sims[j]>=cutoff: assign[i]=j; continue
        centers.append(i); assign[i]=len(centers)-1
    return assign,centers

assign,centers=leader_cluster(fps,0.4)
meta['cluster']=assign

def mk():
    sp=(ytr==0).sum()/ytr.sum()
    return XGBClassifier(n_estimators=500,learning_rate=0.05,max_depth=6,subsample=0.8,
        colsample_bytree=0.8,scale_pos_weight=sp,eval_metric='logloss',random_state=SEED,n_jobs=-1)

rng=np.random.default_rng(42)

def boot_ci(y,yp,pr,n=1000):
    out={'MCC':[], 'BA':[], 'AUC':[]}
    idx=np.arange(len(y))
    for _ in range(n):
        b=rng.choice(idx,len(idx),replace=True)
        if len(np.unique(y[b]))<2: continue
        out['MCC'].append(matthews_corrcoef(y[b],yp[b])); out['BA'].append(balanced_accuracy_score(y[b],yp[b]))
        out['AUC'].append(roc_auc_score(y[b],pr[b]))
    return {k:(float(np.percentile(v,2.5)),float(np.percentile(v,97.5))) for k,v in out.items() if v}

def run(tr_mask,te_mask,name,extra=None):
    m=mk(); m.fit(Xall[tr_mask],yall[tr_mask])
    yp=m.predict(Xall[te_mask]); pr=m.predict_proba(Xall[te_mask])[:,1]; y=yall[te_mask]
    row=dict(protocol=name,n_train=int(tr_mask.sum()),n_test=int(te_mask.sum()),
             test_active_rate=float(y.mean()),
             MCC=matthews_corrcoef(y,yp),BA=balanced_accuracy_score(y,yp),
             AUC=roc_auc_score(y,pr) if len(np.unique(y))>1 else np.nan,
             majority_BA=0.5)
    ci=boot_ci(y,yp,pr)
    for k,(lo,hi) in ci.items(): row[k+'_lo'],row[k+'_hi']=lo,hi
    if extra: row.update(extra)
    return row

rows=[]
rows.append(run((meta.orig=='train').values,(meta.orig=='test').values,'Scaffold split'))
print("scaffold done")
cl=meta.groupby('cluster').size().sort_values(ascending=False)
order=rng.permutation(cl.index.values); cum=0; test_cl=set()
for c in order:
    if cum>=0.20*len(meta): break
    test_cl.add(c); cum+=cl[c]
te_mask=meta.cluster.isin(test_cl).values
rows.append(run(~te_mask,te_mask,'Cluster split (Tanimoto 0.4)'))
print("cluster done: test n =",te_mask.sum())
gb=pd.DataFrame(rows)
print(gb[['protocol','n_test','test_active_rate','MCC','MCC_lo','MCC_hi','BA','AUC']].round(4).to_string(index=False))
gb.to_csv("_gb_partial.csv",index=False)

cts=meta.chemotype.value_counts()
for c in cts.index:
    sel=(meta.chemotype==c).values
    if sel.sum()<30: continue
    if len(np.unique(yall[sel]))<2: continue
    r=run(~sel,sel,f'Leave-out: {c}',extra=dict(held_out_chemotype=c))
    rows.append(r)
    print(f"{c:18s} n={sel.sum():5d} act={yall[sel].mean():.2f} MCC={r['MCC']:.3f} "
          f"[{r.get('MCC_lo',float('nan')):.3f},{r.get('MCC_hi',float('nan')):.3f}] BA={r['BA']:.3f} AUC={r['AUC']:.3f}",flush=True)
gb=pd.DataFrame(rows); gb.to_csv("_gb_partial.csv",index=False)
print("\nSAVED", gb.shape)

rows2=[]; preds={}
def run2(tr_mask,te_mask,name):
    m=mk(); m.fit(Xall[tr_mask],yall[tr_mask])
    yp=m.predict(Xall[te_mask]); pr=m.predict_proba(Xall[te_mask])[:,1]; y=yall[te_mask]
    tn,fp,fn,tp=confusion_matrix(y,yp,labels=[0,1]).ravel()
    preds[name]=dict(y=y,yp=yp,pr=pr,idx=np.where(te_mask)[0])
    return dict(protocol=name,n_test=int(te_mask.sum()),active_rate=float(y.mean()),
        MCC=matthews_corrcoef(y,yp),BA=balanced_accuracy_score(y,yp),
        AUC=roc_auc_score(y,pr) if len(np.unique(y))>1 else np.nan,
        Sens=tp/(tp+fn) if (tp+fn) else np.nan, Spec=tn/(tn+fp) if (tn+fp) else np.nan,
        n_active=int(y.sum()), n_active_found=int(tp))
rows2.append(run2((meta.orig=='train').values,(meta.orig=='test').values,'Scaffold split'))
print(rows2[-1],flush=True)
cl=meta.groupby('cluster').size(); rng2=np.random.default_rng(42)
order=rng2.permutation(cl.index.values); cum=0; test_cl=set()
for c in order:
    if cum>=0.20*len(meta): break
    test_cl.add(c); cum+=cl[c]
tm=meta.cluster.isin(test_cl).values
rows2.append(run2(~tm,tm,'Cluster split (Tanimoto 0.4)')); print(rows2[-1],flush=True)
for c in meta.chemotype.value_counts().index:
    sel=(meta.chemotype==c).values
    if sel.sum()<30 or len(np.unique(yall[sel]))<2: continue
    rows2.append(run2(~sel,sel,f'Leave-out: {c}'))
    r=rows2[-1]; print(f"{c:16s} n={r['n_test']:5d} MCC={r['MCC']:.3f} Sens={r['Sens']:.3f} Spec={r['Spec']:.3f} "
                      f"actives found {r['n_active_found']}/{r['n_active']}",flush=True)
gs=pd.DataFrame(rows2); gs.to_csv("generalization_sens.csv",index=False)

rep_rows=[]
cl2=meta.groupby('cluster').size()
for s in range(5):
    r=np.random.default_rng(100+s); order=r.permutation(cl2.index.values); cum=0; tc=set()
    for c in order:
        if cum>=0.20*len(meta): break
        tc.add(c); cum+=cl2[c]
    tm2=meta.cluster.isin(tc).values
    m=mk(); m.fit(Xall[~tm2],yall[~tm2])
    yp=m.predict(Xall[tm2]); pr=m.predict_proba(Xall[tm2])[:,1]; y=yall[tm2]
    tn,fp,fn,tp=confusion_matrix(y,yp,labels=[0,1]).ravel()
    rep_rows.append(dict(seed=100+s,n_test=int(tm2.sum()),active_rate=float(y.mean()),
        MCC=matthews_corrcoef(y,yp),BA=balanced_accuracy_score(y,yp),AUC=roc_auc_score(y,pr),
        Sens=tp/(tp+fn),Spec=tn/(tn+fp)))
    print(rep_rows[-1],flush=True)
rep=pd.DataFrame(rep_rows); rep.to_csv("cluster_split_repeats.csv",index=False)
print("\nMCC mean %.3f  sd %.3f  range %.3f-%.3f"%(rep.MCC.mean(),rep.MCC.std(),rep.MCC.min(),rep.MCC.max()))

gb2=pd.read_csv("_gb_partial.csv")
m2=gs.merge(gb2[['protocol','MCC_lo','MCC_hi','AUC_lo','AUC_hi']],on='protocol',how='left')
cl_row=dict(protocol='Cluster split (Tanimoto 0.4)',n_test=int(rep.n_test.mean()),active_rate=rep.active_rate.mean(),
        MCC=rep.MCC.mean(),BA=rep.BA.mean(),AUC=rep.AUC.mean(),Sens=rep.Sens.mean(),Spec=rep.Spec.mean(),
        n_active=np.nan,n_active_found=np.nan,MCC_lo=rep.MCC.min(),MCC_hi=rep.MCC.max(),
        AUC_lo=rep.AUC.min(),AUC_hi=rep.AUC.max())
m2=m2[m2.protocol!='Cluster split (Tanimoto 0.4)']
m2=pd.concat([m2,pd.DataFrame([cl_row])],ignore_index=True)
m2['kind']=np.where(m2.protocol.str.startswith('Leave-out'),'Leave-one-chemotype-out','Split protocol')
m2['label']=m2.protocol.str.replace('Leave-out: ','',regex=False)
m2=m2.sort_values(['kind','MCC'],ascending=[True,False])
m2.to_csv("generalization_benchmark.csv",index=False)
print(m2[['protocol','n_test','active_rate','MCC','MCC_lo','MCC_hi','AUC','Sens','Spec','n_active_found','n_active']].round(3).to_string(index=False))