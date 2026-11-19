# AUDIT RECORD for `original_template_enumeration.csv` - automated reconstruction of the producing cell by the platform lineage tool.
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

import os, sys, pickle, warnings
import numpy as np, pandas as pd
warnings.filterwarnings('ignore')
from rdkit import Chem, DataStructs, RDLogger
from rdkit.Chem import Descriptors, QED, rdMolDescriptors, AllChem, Crippen, RDConfig
from rdkit.Chem.Scaffolds import MurckoScaffold
from rdkit.ML.Descriptors import MoleculeDescriptors
from rdkit.DataStructs import ConvertToNumpyArray
RDLogger.DisableLog('rdApp.*')
sys.path.append(os.path.join(RDConfig.RDContribDir, 'SA_Score'))
import sascorer

R = str(DATA_DIR)

ch = pd.read_csv(P("chemotype_assignments.csv"))
CP = pickle.load(open(P("conformal_predictor.pkl"), 'rb'))
MDL, CAL, KEEP, DCOLS = CP['model'], CP['cal'], CP['keep'], CP['desc_cols']

tr = pd.read_csv(R + r"\train.csv")
te = pd.read_csv(R + r"\test.csv")
meta = pd.concat([tr.assign(orig='train'), te.assign(orig='test')], ignore_index=True)
meta = meta.merge(ch[['std_smiles', 'chemotype']].drop_duplicates('std_smiles'), on='std_smiles', how='left')
act = ch[ch.label == 1].copy()

AZOLE = [Chem.MolFromSmarts(s) for s in ["[CX4]n1cncn1", "[CX4]n1ccnc1", "[CX4]n1nccn1", "c1ncnn1", "n1cncn1"]]
def has_azole(m): return any(m.HasSubstructMatch(q) for q in AZOLE)

PAT = [("Triazole", ["[CX4]n1cncn1", "[CX4]n1nccn1"]), ("Imidazole", ["[CX4]n1ccnc1"]),
       ("Fluoropyrimidine", ["Nc1nc(=O)[nH]cc1F", "Fc1c[nX2,nX3]c(=O)[nX2,nX3]c1"]),
       ("Allylamine", ["[NX3;!$(NC=O)](C)CC=C"]), ("Thiocarbamate", ["[NX3]C(=S)O", "[NX3]C(=O)S"]),
       ("Morpholine", ["C1COCCN1"]), ("Benzimidazole", ["c1ccc2[nH]cnc2c1"])]
COMP = [(n, [Chem.MolFromSmarts(s) for s in p]) for n, p in PAT]

def chemotype(smi):
    m = Chem.MolFromSmiles(smi)
    if m is None: return "Unparsed"
    mw = Descriptors.MolWt(m)
    rings = [len(r) for r in m.GetRingInfo().AtomRings()]
    if mw > 800 and rdMolDescriptors.CalcNumAmideBonds(m) >= 5: return "Echinocandin"
    if rings and max(rings) >= 20 and mw > 600: return "Polyene"
    for n, ms in COMP:
        if any(m.HasSubstructMatch(q) for q in ms): return n
    return "Other"

GEN = AllChem.GetMorganGenerator(radius=2, fpSize=2048)
TRAIN_FPS = [GEN.GetFingerprint(m) for m in (Chem.MolFromSmiles(s) for s in meta.std_smiles) if m]
TRAIN_SCAF = set(meta.std_smiles.apply(lambda s: MurckoScaffold.MurckoScaffoldSmiles(smiles=s, includeChirality=False)))

ALL_SUBS = ['F', 'Cl', 'Br', 'C', 'CC', 'CCC', 'C(C)C', 'OC', 'OCC', 'OC(C)C', 'O', 'C(F)(F)F', 'C#N', 'S(=O)(=O)N',
            'C(=O)N', 'C(=O)OC', 'N(=O)=O', 'n1ccnc1', 'n1cccc1', 'c1ccncc1', 'c1cncnc1', 'c1ccoc1', 'c1ccsc1',
            'N', 'NC', 'NCC', 'NC(C)C', 'NC(=O)C']
TEMPLATES = ['OC(Cn1cncn1)(Cn1cncn1)c1ccc({R1})cc1', 'OC(Cn1cncn1)(Cn1cncn1)c1cc({R1})cc(F)c1',
             'OC(Cn1cncn1)(Cn1cncn1)c1ccc(CC{R1})cc1', 'OC(Cn1cncn1)(Cc1ncc({R1})cn1)c1ccc(F)cc1',
             'OC(Cn1cncn1)(Cn1cncn1)c1cc(F)c({R1})cc1F']
orig = set()
for t in TEMPLATES:
    for s in ALL_SUBS:
        m = Chem.MolFromSmiles(t.replace('{R1}', s))
        if m: orig.add(Chem.MolToSmiles(m))
orig = sorted(orig)
rows = []
for s in orig:
    m = Chem.MolFromSmiles(s)
    fp = GEN.GetFingerprint(m)
    mx = float(max(DataStructs.BulkTanimotoSimilarity(fp, TRAIN_FPS)))
    scaf = MurckoScaffold.MurckoScaffoldSmiles(mol=m, includeChirality=False)
    rows.append(dict(smiles=s, max_tanimoto=mx, novelty=1-mx, sa=sascorer.calculateScore(m),
                     qed=QED.qed(m), scaffold=scaf, scaffold_novel=int(scaf not in TRAIN_SCAF),
                     chemotype=chemotype(s)))
od = pd.DataFrame(rows)
od.to_csv("original_template_enumeration.csv", index=False)