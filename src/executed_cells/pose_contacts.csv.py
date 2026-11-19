# AUDIT RECORD for `pose_contacts.csv` - raw producing cell (automated reconstruction failed validation and the platform returned the raw cell).
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

# [lineage] reconstruction failed validation (max_tokens) — raw producing cell below; see Execution Log for full trace

import os, numpy as np, pandas as pd
D="dock"
cs=pd.read_csv(os.path.join(D,"candidate_scores.csv")).dropna(subset=['affinity'])
cand=pd.read_csv("docking_candidates.csv")
m=cand.merge(cs,on='cand_id',how='left')
ctrl=cs[cs.cand_id.isin(['fluconazole','posaconazole'])]
print("controls:"); print(ctrl.to_string(index=False))
print("\ncandidates docked:",m.affinity.notna().sum(),"| median %.2f | best %.2f"%(m.affinity.median(),m.affinity.min()))
m.to_csv("candidates_with_docking.csv",index=False)

def pose_atoms(path):
    """first MODEL heavy atoms from a vina pdbqt"""
    out=[]; started=False
    for l in open(path):
        if l.startswith("MODEL"):
            if started: break
            started=True; continue
        if l.startswith("ENDMDL"): break
        if l.startswith(("ATOM","HETATM")):
            if l[12:16].strip().startswith("H"): continue
            out.append(l)
    return out
rec_lines=[l for l in open(os.path.join(D,"receptor_A.pdb")) if l.startswith(("ATOM","HETATM"))]
rec_xyz=np.array([[float(l[30:38]),float(l[38:46]),float(l[46:54])] for l in rec_lines])
rec_res=[(l[17:20].strip(),l[22:26].strip()) for l in rec_lines]
fe=[l for l in rec_lines if l[17:20]=="HEM" and l[12:16].strip()=="FE"][0]
fexyz=np.array([float(fe[30:38]),float(fe[38:46]),float(fe[46:54])])

def build_complex(tag,src,dest):
    pa=pose_atoms(src)
    lig=[]
    for i,l in enumerate(pa):
        el=(l[76:78].strip() or l[12:16].strip()[0])
        lig.append("HETATM%5d %-4s LIG L 901    %8.3f%8.3f%8.3f  1.00  0.00          %2s\n"%(
            9000+i, l[12:16].strip()[:4], float(l[30:38]),float(l[38:46]),float(l[46:54]), el[:2]))
    with open(dest,"w") as f:
        f.write("REMARK  docked complex: C. albicans CYP51 (5TZ1 chain A + heme) with %s\n"%tag)
        f.writelines(rec_lines); f.write("TER\n"); f.writelines(lig); f.write("TER\nEND\n")
    xyz=np.array([[float(l[30:38]),float(l[38:46]),float(l[46:54])] for l in pa])
    d=np.linalg.norm(rec_xyz[:,None,:]-xyz[None,:,:],axis=2)
    near=sorted({f"{rec_res[i][0]}{rec_res[i][1]}" for i in np.where(d.min(1)<4.0)[0]})
    dfe=float(np.linalg.norm(xyz-fexyz,axis=1).min())
    return near,dfe,len(pa)

top=m.dropna(subset=['affinity']).sort_values('affinity').head(5)
rows=[]
os.makedirs("complexes",exist_ok=True)
for _,r in top.iterrows():
    src=os.path.join(D,"cand_poses",f"{r.cand_id}_out.pdbqt")
    dest=f"complexes/{r.cand_id}_CYP51_complex.pdb"
    near,dfe,na=build_complex(r.cand_id,src,dest)
    rows.append(dict(cand_id=r.cand_id,affinity=r.affinity,n_heavy=na,min_dist_to_heme_Fe=round(dfe,2),
                     n_contacts=len(near),contacts="; ".join(near)))
    print(f"{r.cand_id} {r.affinity:7.2f}  Fe {dfe:5.2f} A  contacts {len(near)}")
for c in ['fluconazole','posaconazole']:
    src=os.path.join(D,"cand_poses",f"{c}_out.pdbqt")
    if os.path.exists(src):
        near,dfe,na=build_complex(c,src,f"complexes/{c}_CYP51_complex.pdb")
        a=float(cs[cs.cand_id==c].affinity.iloc[0])
        rows.append(dict(cand_id=c,affinity=a,n_heavy=na,min_dist_to_heme_Fe=round(dfe,2),
                         n_contacts=len(near),contacts="; ".join(near)))
        print(f"{c:14s} {a:7.2f}  Fe {dfe:5.2f} A  contacts {len(near)}")
pc=pd.DataFrame(rows); pc.to_csv("pose_contacts.csv",index=False)
print("\ncomplex files:",sorted(os.listdir("complexes")))
