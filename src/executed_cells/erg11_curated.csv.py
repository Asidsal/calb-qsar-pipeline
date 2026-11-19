# AUDIT RECORD for `erg11_curated.csv` - automated reconstruction of the producing cell by the platform lineage tool.
# Absolute local paths were replaced by P(<filename>); nothing else was edited. This file has NOT been re-run.
# Inputs recorded: chembl_erg11_raw.pkl, 01_data_curation.ipynb
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
    # unit harmonisation -> nM
    def to_nM(r):
        u=r.standard_units; v=r.standard_value
        if u=='nM': return v
        if u=='uM': return v*1e3
        if u in ('ug.mL-1','ug ml-1'):
            m=Chem.MolFromSmiles(r.canonical_smiles)
            if m is None: return np.nan
            return v/Descriptors.MolWt(m)*1e6   # ug/mL -> uM -> nM
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

PATTERNS=[
 ("Triazole",        ["[CX4]n1cncn1","[CX4]n1nccn1"]),
 ("Imidazole",       ["[CX4]n1ccnc1"]),
 ("Fluoropyrimidine",["Nc1nc(=O)[nH]cc1F","Fc1c[nX2,nX3]c(=O)[nX2,nX3]c1"]),
 ("Allylamine",      ["[NX3;!$(NC=O)](C)CC=C"]),
 ("Thiocarbamate",   ["[NX3]C(=S)O","[NX3]C(=O)S"]),
 ("Morpholine",      ["C1COCCN1"]),
 ("Benzimidazole",   ["c1ccc2[nH]cnc2c1"]),
]
COMPILED=[]
for n,p in PATTERNS:
    ms=[Chem.MolFromSmarts(s) for s in p]
    bad=[s for s,m in zip(p,ms) if m is None]
    assert not bad, (n,bad)
    COMPILED.append((n,ms))

def chemotype(smi):
    m=Chem.MolFromSmiles(smi)
    if m is None: return "Unparsed"
    mw=Descriptors.MolWt(m); ri=m.GetRingInfo()
    rings=[len(r) for r in ri.AtomRings()]
    from rdkit.Chem import rdMolDescriptors
    if mw>800 and rdMolDescriptors.CalcNumAmideBonds(m)>=5: return "Echinocandin"
    if rings and max(rings)>=20 and mw>600: return "Polyene"
    for n,ms in COMPILED:
        if any(m.HasSubstructMatch(q) for q in ms): return n
    return "Other"

erg['chemotype']=erg.std_smiles.apply(chemotype)
erg.to_csv("erg11_curated.csv",index=False)