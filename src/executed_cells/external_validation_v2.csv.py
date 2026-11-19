# AUDIT RECORD for `external_validation_v2.csv` - automated reconstruction of the producing cell by the platform lineage tool.
# Absolute local paths were replaced by P(<filename>); nothing else was edited. This file has NOT been re-run.
# Inputs recorded: conformal_predictor.pkl, chemotype_assignments.csv, X_train_desc.csv
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

import pandas as pd
import numpy as np
import pickle
import requests
import time
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem.MolStandardize import rdMolStandardize
from rdkit.Chem import Descriptors, rdMolDescriptors
from rdkit.ML.Descriptors import MoleculeDescriptors
from rdkit.DataStructs import ConvertToNumpyArray

# Load conformal predictor
CP = pickle.load(open(P("conformal_predictor.pkl"), 'rb'))
desc_cols = CP['desc_cols']
mdl = CP['model']
cal = CP['cal']
keep = CP['keep']

calc = MoleculeDescriptors.MolecularDescriptorCalculator(desc_cols)
med = pd.read_csv(P("X_train_desc.csv")).median()
gen = AllChem.GetMorganGenerator(radius=2, fpSize=2048)

# Load chemotype assignments for pool membership check
ch = pd.read_csv(P("chemotype_assignments.csv"))

# SMARTS patterns for chemotype assignment
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
    COMPILED.append((n, ms))


def chemotype(smi):
    m = Chem.MolFromSmiles(smi)
    if m is None:
        return "Unparsed"
    mw = Descriptors.MolWt(m)
    ri = m.GetRingInfo()
    rings = [len(r) for r in ri.AtomRings()]
    if mw > 800 and rdMolDescriptors.CalcNumAmideBonds(m) >= 5:
        return "Echinocandin"
    if rings and max(rings) >= 20 and mw > 600:
        return "Polyene"
    for n, ms in COMPILED:
        if any(m.HasSubstructMatch(q) for q in ms):
            return n
    return "Other"


def std(s):
    try:
        m = Chem.MolFromSmiles(s)
        if m is None:
            return None
        m = rdMolStandardize.LargestFragmentChooser().choose(m)
        m = rdMolStandardize.Uncharger().uncharge(m)
        return Chem.MolToSmiles(m, canonical=True)
    except:
        return None


def ikey(s):
    m = Chem.MolFromSmiles(s)
    return Chem.MolToInchiKey(m).split('-')[0] if m else None


# Fetch compounds from PubChem
S = requests.Session()
PUG = "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/{}/property/SMILES,ConnectivitySMILES,MolecularFormula/JSON"


def lookup(name):
    for _ in range(3):
        try:
            r = S.get(PUG.format(requests.utils.quote(name)), timeout=30)
            if r.status_code == 404:
                return None
            r.raise_for_status()
            p = r.json()["PropertyTable"]["Properties"][0]
            smi = p.get("SMILES") or p.get("ConnectivitySMILES")
            if smi:
                return dict(query=name, cid=p["CID"], smiles=smi, formula=p.get("MolecularFormula"))
        except Exception:
            time.sleep(1.0)
    return None


ANTIFUNGAL = {
    'Triazole': ['fluconazole', 'itraconazole', 'voriconazole', 'posaconazole', 'isavuconazole', 'ravuconazole',
                 'albaconazole', 'efinaconazole', 'terconazole', 'propiconazole', 'oteseconazole', 'fosravuconazole'],
    'Imidazole': ['ketoconazole', 'miconazole', 'econazole', 'clotrimazole', 'butoconazole', 'sulconazole',
                  'oxiconazole', 'sertaconazole', 'tioconazole', 'bifonazole', 'isoconazole', 'fenticonazole', 'luliconazole'],
    'Polyene': ['amphotericin B', 'nystatin', 'natamycin'],
    'Echinocandin': ['caspofungin', 'micafungin', 'anidulafungin', 'rezafungin'],
    'Allylamine': ['terbinafine', 'naftifine', 'butenafine'],
    'Other antifungal': ['flucytosine', 'griseofulvin', 'ciclopirox', 'amorolfine', 'tolnaftate', 'tavaborole',
                         'ibrexafungerp', 'olorofim', 'manogepix', 'haloprogin', 'undecylenic acid'],
}
DECOY = {
    'Statin': ['atorvastatin', 'simvastatin', 'lovastatin', 'fluvastatin', 'pravastatin', 'rosuvastatin'],
    'Cardiovascular': ['amlodipine', 'verapamil', 'lisinopril', 'losartan', 'atenolol', 'propranolol'],
    'Antidepressant': ['sertraline', 'fluoxetine', 'amitriptyline', 'citalopram', 'venlafaxine', 'bupropion'],
    'Antipsychotic': ['chlorpromazine', 'thioridazine', 'haloperidol', 'olanzapine', 'risperidone', 'quetiapine'],
    'NSAID/analgesic': ['ibuprofen', 'diclofenac', 'naproxen', 'celecoxib', 'acetaminophen', 'aspirin'],
    'Antibacterial': ['ciprofloxacin', 'amoxicillin', 'doxycycline', 'azithromycin', 'trimethoprim', 'linezolid'],
    'Gastrointestinal': ['omeprazole', 'lansoprazole', 'ranitidine', 'famotidine', 'ondansetron', 'metoclopramide'],
    'Antihistamine': ['diphenhydramine', 'promethazine', 'loratadine', 'cetirizine', 'fexofenadine', 'chlorpheniramine'],
    'Antidiabetic': ['metformin', 'glibenclamide', 'pioglitazone', 'sitagliptin', 'gliclazide', 'repaglinide'],
    'Other non-antifungal': ['ivermectin', 'albendazole', 'warfarin', 'furosemide', 'levothyroxine', 'prednisolone'],
}

rows = []
miss = []
for panel, D in [('Antifungal', ANTIFUNGAL), ('Non-antifungal', DECOY)]:
    for grp, names in D.items():
        for nm in names:
            r = lookup(nm)
            if r:
                r.update(panel=panel, group=grp)
                rows.append(r)
            else:
                miss.append((panel, nm))

p = pd.DataFrame(rows)
p['std_smiles'] = p.smiles.apply(std)
p = p[p.std_smiles.notna()].drop_duplicates('std_smiles').reset_index(drop=True)

# Check training pool membership using InChIKey
pool_ikeys = set(ch.std_smiles.apply(ikey).dropna())

p['ikey'] = p.std_smiles.apply(ikey)

def place(k):
    return 'TRAIN' if k in pool_ikeys else 'external'

p['placement'] = p.ikey.apply(place)
p = p.drop_duplicates('ikey').reset_index(drop=True)
p['chemotype'] = p.std_smiles.apply(chemotype)

# Assign expected labels
CANDIDA_ACTIVE_CLASSES = {'Triazole', 'Imidazole', 'Polyene', 'Echinocandin'}
NON_CANDIDA_ANTIFUNGAL = {'terbinafine', 'tolnaftate', 'olorofim'}
UNSCORED = {'tavaborole', 'undecylenic acid'}


def lab(r):
    if r['query'] in UNSCORED:
        return np.nan, 'equivocal anti-Candida spectrum - not scored'
    if r.panel == 'Non-antifungal':
        return 0, 'no antifungal indication'
    if r['query'] in NON_CANDIDA_ANTIFUNGAL:
        return 0, 'antifungal but not anti-Candida spectrum'
    if r.group in CANDIDA_ACTIVE_CLASSES:
        return 1, f'{r.group} - established anti-Candida class'
    return np.nan, 'antifungal of uncertain anti-Candida spectrum - not scored'


p[['expected', 'label_basis']] = p.apply(lambda r: pd.Series(lab(r)), axis=1)

# Filter to external, scored compounds
ext = p[(p.placement == 'external') & (p.expected.notna())].reset_index(drop=True)


def featurise(smiles):
    F = []
    D = []
    for s in smiles:
        m = Chem.MolFromSmiles(s)
        arr = np.zeros(2048, dtype=np.uint8)
        if m is not None:
            ConvertToNumpyArray(gen.GetFingerprint(m), arr)
            d = list(calc.CalcDescriptors(m))
            d = [med.iloc[i] if (v != v or abs(v) > 1e10) else v for i, v in enumerate(d)]
        else:
            d = list(med.values)
        F.append(arr)
        D.append(d)
    return np.hstack([np.array(F), np.array(D)])


def pvals(P, groups, cal):
    pv = np.ones((len(P), 2))
    for l in (0, 1):
        a = 1.0 - P[:, l]
        for g in np.unique(groups):
            mk_ = groups == g
            ref = cal.get((g, l))
            if ref is None or len(ref) == 0:
                continue
            pv[mk_, l] = (len(ref) - np.searchsorted(ref, a[mk_], side='left') + 1) / (len(ref) + 1)
    return pv


X = featurise(ext.std_smiles.tolist())
P = mdl.predict_proba(X)
g = np.array([c if c in keep else 'POOLED' for c in ext.chemotype])
pv = pvals(P, g, cal)

ext['prob_active'] = P[:, 1]
ext['p_inactive'] = pv[:, 0]
ext['p_active'] = pv[:, 1]
ext['point_pred'] = (P[:, 1] >= 0.5).astype(int)
ext['point_correct'] = (ext.point_pred == ext.expected).astype(int)

for e in (0.01, 0.05, 0.10, 0.20):
    reg = pv > e
    single = reg.sum(1) == 1
    lab_arr = np.where(single, reg.argmax(1), -1)
    ext[f'call_eps{e}'] = lab_arr

ext.to_csv('external_validation_v2.csv', index=False)