# AUDIT RECORD for `chemotype_assignments.csv` - automated reconstruction of the producing cell by the platform lineage tool.
# Absolute local paths were replaced by P(<filename>); nothing else was edited. This file has NOT been re-run.
# Inputs recorded: CAlbicans_curated.csv
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
from rdkit.Chem import Descriptors, rdMolDescriptors

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
    if mw>800 and rdMolDescriptors.CalcNumAmideBonds(m)>=5: return "Echinocandin"
    if rings and max(rings)>=20 and mw>600: return "Polyene"
    for n,ms in COMPILED:
        if any(m.HasSubstructMatch(q) for q in ms): return n
    return "Other"

phen=pd.read_csv(P("CAlbicans_curated.csv"))
phen['chemotype']=phen.std_smiles.apply(chemotype)
phen.to_csv("chemotype_assignments.csv",index=False)