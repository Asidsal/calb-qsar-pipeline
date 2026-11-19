# AUDIT RECORD for `candidates_scored.csv` - raw producing cell (automated reconstruction failed validation and the platform returned the raw cell).
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

# [lineage] reconstruction failed validation (max_tokens) — raw producing cell below; see Execution Log for full trace

import pandas as pd, numpy as np
ds=pd.read_csv("docking_scores_enrichment.csv")
print("precision of a docking cut-off, measured on the enrichment set:")
rows=[]
for t in [-9.5,-10.0,-10.5,-11.0,-11.5]:
    sel=ds[ds.affinity<=t]
    if len(sel)<5: continue
    rows.append(dict(threshold=t,n=len(sel),actives=int(sel.y.sum()),
                     precision=sel.y.mean(),recall=sel.y.sum()/ds.y.sum()))
pr=pd.DataFrame(rows); pr['baseline']=ds.y.mean()
print(pr.round(3).to_string(index=False))
GATE=-10.5
print(f"\nchosen plausibility gate: affinity <= {GATE} kcal/mol "
      f"(precision {float(pr[pr.threshold==GATE].precision.iloc[0]):.2f} vs baseline {ds.y.mean():.2f})")
c=pd.read_csv("candidates_with_docking.csv")
h=pd.read_csv("dock/human_scores.csv").rename(columns={'cand_id':'cand_id'}).dropna(subset=['affinity_human'])
m=c.merge(h[['cand_id','affinity_human']],on='cand_id',how='left')
m['selectivity']=m.affinity-m.affinity_human        # negative = favours fungal enzyme
m=m.dropna(subset=['affinity','affinity_human'])
print("\ncandidates with both scores:",len(m))
print("selectivity (fungal - human, kcal/mol): median %.2f | favouring fungal %d | favouring human %d"%(
    m.selectivity.median(),(m.selectivity<0).sum(),(m.selectivity>0).sum()))
m['docking_supported']=(m.affinity<=GATE).astype(int)
m['selectivity_favourable']=(m.selectivity<=-0.5).astype(int)
m['evidence']=np.where(m.docking_supported==1,'QSAR + docking','QSAR only')
print("\nevidence classes:\n",m.evidence.value_counts().to_string())
print("\nQSAR+docking AND selectivity-favourable:",int(((m.docking_supported==1)&(m.selectivity_favourable==1)).sum()))
m.to_csv("candidates_scored.csv",index=False)
print("\ntop by combined evidence:")
t=m[(m.docking_supported==1)].sort_values('selectivity')
print(t[['cand_id','track','prob_active','affinity','affinity_human','selectivity','max_tanimoto','sa','qed']].head(12).round(3).to_string(index=False))
