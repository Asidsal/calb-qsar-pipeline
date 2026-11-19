# AUDIT RECORD for `dual_axis_overlap.csv` - automated reconstruction of the producing cell by the platform lineage tool.
# Absolute local paths were replaced by P(<filename>); nothing else was edited. This file has NOT been re-run.
# Inputs recorded: chembl_erg11_raw.pkl, CAlbicans_curated.csv
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

import pandas as pd, numpy as np
from rdkit import Chem
from rdkit.Chem import Descriptors
from rdkit.Chem.MolStandardize import rdMolStandardize
from rdkit.Chem.Scaffolds import MurckoScaffold
from scipy.stats import gmean

def standardize(smi):
    try:
        m=Chem.MolFromSmiles(smi)
        if m is None: return None
        m=rdMolStandardize.LargestFragmentChooser().choose(m)
        m=rdMolStandardize.Uncharger().uncharge(m)
        return Chem.MolToSmiles(m,canonical=True)
    except: return None

def curate_target(df,name):
    d=df[df.standard_type.isin(['IC50','Ki','Kd','EC50'])].copy()
    d=d[d.standard_relation=='=']
    d=d[d.canonical_smiles.notna() & d.standard_value.notna()]
    d['standard_value']=pd.to_numeric(d.standard_value,errors='coerce')
    def to_nM(r):
        u=r.standard_units; v=r.standard_value
        if u=='nM': return v
        if u=='uM': return v*1e3
        if u in ('ug.mL-1','ug ml-1'):
            m=Chem.MolFromSmiles(r.canonical_smiles)
            if m is None: return np.nan
            return v/Descriptors.MolWt(m)*1e6
        return np.nan
    d['nM']=d.apply(to_nM,axis=1)
    d=d[d.nM.notna() & (d.nM>0)]
    d['std_smiles']=d.canonical_smiles.apply(standardize)
    d=d[d.std_smiles.notna()]
    g=d.groupby('molecule_chembl_id').agg(std_smiles=('std_smiles','first'),
        nM=('nM',lambda x: float(gmean(x))), n_records=('nM','count'),
        types=('standard_type',lambda x:'|'.join(sorted(set(x))))).reset_index()
    g['pActivity']=-np.log10(g.nM*1e-9)
    g['scaffold']=g.std_smiles.apply(lambda s: MurckoScaffold.MurckoScaffoldSmiles(smiles=s,includeChirality=False))
    print(f"{name}: {len(g)} compounds | pActivity {g.pActivity.min():.2f}-{g.pActivity.max():.2f} median {g.pActivity.median():.2f}")
    return g

ca=pd.read_pickle(P("chembl_erg11_raw.pkl"))
erg=curate_target(ca,"Erg11 (C. albicans)")

phen=pd.read_csv(P("CAlbicans_curated.csv"))
ov=phen.merge(erg[['std_smiles','pActivity','nM']],on='std_smiles',how='inner')
print(f"\nOverlap phenotypic (n={len(phen)}) x Erg11 (n={len(erg)}): {len(ov)} compounds")
ov.to_csv("dual_axis_overlap.csv",index=False)