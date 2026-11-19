# AUDIT RECORD for `external_panel_v2.csv` - automated reconstruction of the producing cell by the platform lineage tool.
# Absolute local paths were replaced by P(<filename>); nothing else was edited. This file has NOT been re-run.
# Inputs recorded: conformal_predictor.pkl, partA_antifungal_validation.csv, train.csv, test.csv
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

import requests, pandas as pd, time
from rdkit import Chem
from rdkit.Chem.MolStandardize import rdMolStandardize

S=requests.Session()
PUG="https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/{}/property/SMILES,ConnectivitySMILES,MolecularFormula/JSON"
def lookup(name):
    for _ in range(3):
        try:
            r=S.get(PUG.format(requests.utils.quote(name)),timeout=30)
            if r.status_code==404: return None
            r.raise_for_status()
            p=r.json()["PropertyTable"]["Properties"][0]
            smi=p.get("SMILES") or p.get("ConnectivitySMILES")
            if smi: return dict(query=name,cid=p["CID"],smiles=smi,formula=p.get("MolecularFormula"))
        except Exception: time.sleep(1.0)
    return None

ANTIFUNGAL={
 'Triazole':['fluconazole','itraconazole','voriconazole','posaconazole','isavuconazole','ravuconazole',
             'albaconazole','efinaconazole','terconazole','propiconazole','oteseconazole','fosravuconazole'],
 'Imidazole':['ketoconazole','miconazole','econazole','clotrimazole','butoconazole','sulconazole',
              'oxiconazole','sertaconazole','tioconazole','bifonazole','isoconazole','fenticonazole','luliconazole'],
 'Polyene':['amphotericin B','nystatin','natamycin'],
 'Echinocandin':['caspofungin','micafungin','anidulafungin','rezafungin'],
 'Allylamine':['terbinafine','naftifine','butenafine'],
 'Other antifungal':['flucytosine','griseofulvin','ciclopirox','amorolfine','tolnaftate','tavaborole',
                     'ibrexafungerp','olorofim','manogepix','haloprogin','undecylenic acid'],
}
DECOY={
 'Statin':['atorvastatin','simvastatin','lovastatin','fluvastatin','pravastatin','rosuvastatin'],
 'Cardiovascular':['amlodipine','verapamil','lisinopril','losartan','atenolol','propranolol'],
 'Antidepressant':['sertraline','fluoxetine','amitriptyline','citalopram','venlafaxine','bupropion'],
 'Antipsychotic':['chlorpromazine','thioridazine','haloperidol','olanzapine','risperidone','quetiapine'],
 'NSAID/analgesic':['ibuprofen','diclofenac','naproxen','celecoxib','acetaminophen','aspirin'],
 'Antibacterial':['ciprofloxacin','amoxicillin','doxycycline','azithromycin','trimethoprim','linezolid'],
 'Gastrointestinal':['omeprazole','lansoprazole','ranitidine','famotidine','ondansetron','metoclopramide'],
 'Antihistamine':['diphenhydramine','promethazine','loratadine','cetirizine','fexofenadine','chlorpheniramine'],
 'Antidiabetic':['metformin','glibenclamide','pioglitazone','sitagliptin','gliclazide','repaglinide'],
 'Other non-antifungal':['ivermectin','albendazole','warfarin','furosemide','levothyroxine','prednisolone'],
}
rows=[]; miss=[]
for panel,D in [('Antifungal',ANTIFUNGAL),('Non-antifungal',DECOY)]:
    for grp,names in D.items():
        for nm in names:
            r=lookup(nm)
            if r: r.update(panel=panel,group=grp); rows.append(r)
            else: miss.append((panel,nm))
p=pd.DataFrame(rows)

def std(s):
    try:
        m=Chem.MolFromSmiles(s)
        if m is None: return None
        m=rdMolStandardize.LargestFragmentChooser().choose(m)
        m=rdMolStandardize.Uncharger().uncharge(m)
        return Chem.MolToSmiles(m,canonical=True)
    except: return None

p['std_smiles']=p.smiles.apply(std)
p=p[p.std_smiles.notna()].drop_duplicates('std_smiles').reset_index(drop=True)

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
for n,pts in PATTERNS:
    ms=[Chem.MolFromSmarts(s) for s in pts]
    COMPILED.append((n,ms))

from rdkit.Chem import Descriptors, rdMolDescriptors

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

tr=pd.read_csv(P("train.csv"))
te=pd.read_csv(P("test.csv"))
meta=pd.concat([tr.assign(orig='train'),te.assign(orig='test')],ignore_index=True)

def ikey(s):
    m=Chem.MolFromSmiles(s)
    return Chem.MolToInchiKey(m).split('-')[0] if m else None

meta['ikey']=meta.std_smiles.apply(ikey)
tr_k=set(meta.loc[meta.orig=='train','ikey'].dropna())
te_k=set(meta.loc[meta.orig=='test','ikey'].dropna())

def place(k): return 'TRAIN' if k in tr_k else ('test' if k in te_k else 'external')

p['ikey']=p.std_smiles.apply(ikey)
p['placement']=p.ikey.apply(place)
p=p.drop_duplicates('ikey').reset_index(drop=True)
p['chemotype']=p.std_smiles.apply(chemotype)

import numpy as np

CANDIDA_ACTIVE_CLASSES={'Triazole','Imidazole','Polyene','Echinocandin'}
NON_CANDIDA_ANTIFUNGAL={'terbinafine','tolnaftate','olorofim'}
UNSCORED={'tavaborole','undecylenic acid'}

def lab(r):
    if r['query'] in UNSCORED: return np.nan,'equivocal anti-Candida spectrum - not scored'
    if r.panel=='Non-antifungal': return 0,'no antifungal indication'
    if r['query'] in NON_CANDIDA_ANTIFUNGAL: return 0,'antifungal but not anti-Candida spectrum'
    if r.group in CANDIDA_ACTIVE_CLASSES: return 1,f'{r.group} - established anti-Candida class'
    return np.nan,'antifungal of uncertain anti-Candida spectrum - not scored'

p[['expected','label_basis']]=p.apply(lambda r: pd.Series(lab(r)),axis=1)

p.to_csv("external_panel_v2.csv",index=False)