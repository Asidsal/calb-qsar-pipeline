# AUDIT RECORD for `chem_space_umap.csv` - automated reconstruction of the producing cell by the platform lineage tool.
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

import numpy as np, pandas as pd, time, warnings, os, sys, pickle, re
warnings.filterwarnings('ignore')
from rdkit import Chem, DataStructs, RDLogger
from rdkit.Chem import BRICS, Descriptors, QED, rdMolDescriptors, AllChem, Crippen, RDConfig
from rdkit.Chem.Scaffolds import MurckoScaffold
from rdkit.Chem.FilterCatalog import FilterCatalog, FilterCatalogParams
from rdkit.ML.Descriptors import MoleculeDescriptors
from rdkit.DataStructs import ConvertToNumpyArray
RDLogger.DisableLog('rdApp.*')
sys.path.append(os.path.join(RDConfig.RDContribDir, 'SA_Score'))
import sascorer
import umap

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
    mw = Descriptors.MolWt(m); rings = [len(r) for r in m.GetRingInfo().AtomRings()]
    if mw > 800 and rdMolDescriptors.CalcNumAmideBonds(m) >= 5: return "Echinocandin"
    if rings and max(rings) >= 20 and mw > 600: return "Polyene"
    for n, ms in COMP:
        if any(m.HasSubstructMatch(q) for q in ms): return n
    return "Other"

pp = FilterCatalogParams()
for c in (FilterCatalogParams.FilterCatalogs.PAINS, FilterCatalogParams.FilterCatalogs.BRENK): pp.AddCatalog(c)
CATALOG = FilterCatalog(pp)
CALC = MoleculeDescriptors.MolecularDescriptorCalculator(DCOLS)
MED = pd.read_csv(R + r"\X_train_desc.csv").median()
GEN = AllChem.GetMorganGenerator(radius=2, fpSize=2048)
TRAIN_FPS = [GEN.GetFingerprint(m) for m in (Chem.MolFromSmiles(s) for s in meta.std_smiles) if m]
TRAIN_SCAF = set(meta.std_smiles.apply(lambda s: MurckoScaffold.MurckoScaffoldSmiles(smiles=s, includeChirality=False)))

import random
random.seed(42)
samp = act.groupby('chemotype', group_keys=False).apply(lambda d: d.sample(min(len(d), 400), random_state=42))
frags = set()
for s in samp.std_smiles:
    m = Chem.MolFromSmiles(s)
    if m is None: continue
    try: frags.update(BRICS.BRICSDecompose(m, minFragmentSize=3))
    except Exception: pass
frags = {f for f in frags if f}

def frag_ok(f, allow_azole):
    m = Chem.MolFromSmiles(f)
    if m is None: return False
    h = m.GetNumHeavyAtoms() - sum(1 for a in m.GetAtoms() if a.GetAtomicNum() == 0)
    return (3 <= h <= 22) and (allow_azole or not has_azole(m))

FRAGS_ALL = sorted(f for f in frags if frag_ok(f, True))
FRAGS_NOAZ = sorted(f for f in frags if frag_ok(f, False))

def to_single_ap(fs):
    out = []
    for f in fs:
        if f.count('*') != 1: continue
        m = Chem.MolFromSmiles(re.sub(r'\[\d+\*\]', '[1*]', f))
        if m is None: continue
        if 2 <= m.GetNumHeavyAtoms() - 1 <= 18: out.append(m)
    return out

SIDE_ALL = to_single_ap(FRAGS_ALL)
SIDE_NOAZ = to_single_ap(FRAGS_NOAZ)
JOIN = AllChem.ReactionFromSmarts('[*:1]-[1*].[1*]-[*:2]>>[*:1]-[*:2]')
CUTSMA = Chem.MolFromSmarts('[*]-;!@[*]')

def cut(mol):
    hits = mol.GetSubstructMatches(CUTSMA)
    if not hits: return None
    for _ in range(6):
        a, b = random.choice(hits)
        bd = mol.GetBondBetweenAtoms(a, b)
        if bd is None: continue
        try:
            fr = Chem.FragmentOnBonds(mol, [bd.GetIdx()], addDummies=True, dummyLabels=[(1, 1)])
            parts = Chem.GetMolFrags(fr, asMols=True, sanitizeFrags=True)
        except Exception: continue
        if len(parts) == 2: return parts
    return None

def join(f1, f2):
    try: prods = JOIN.RunReactants((f1, f2))
    except Exception: return None
    for p in prods:
        m = p[0]
        try:
            Chem.SanitizeMol(m)
            smi = Chem.MolToSmiles(m)
            if '*' not in smi: return smi
        except Exception: continue
    return None

def crossover(m1, m2):
    for _ in range(5):
        c1, c2 = cut(m1), cut(m2)
        if c1 and c2:
            s = join(random.choice(c1), random.choice(c2))
            if s: return s
    return None

def mutate(m, side):
    for _ in range(5):
        c = cut(m)
        if c:
            s = join(random.choice(c), random.choice(side))
            if s: return s
    return None

def props(m):
    return dict(MW=Descriptors.MolWt(m), LogP=Crippen.MolLogP(m), HBD=rdMolDescriptors.CalcNumHBD(m),
                HBA=rdMolDescriptors.CalcNumHBA(m), TPSA=Descriptors.TPSA(m),
                RotB=rdMolDescriptors.CalcNumRotatableBonds(m), Heavy=m.GetNumHeavyAtoms(),
                Rings=rdMolDescriptors.CalcNumRings(m))

_cache = {}
def score(smi, no_azole):
    if smi in _cache: return _cache[smi]
    m = Chem.MolFromSmiles(smi)
    if m is None: _cache[smi] = None; return None
    p = props(m)
    bad = (not 250 <= p['MW'] <= 550) or (not -1.0 <= p['LogP'] <= 5.5) or p['HBD'] > 5 or p['HBA'] > 10 \
        or p['RotB'] > 10 or p['TPSA'] > 140 or (not 17 <= p['Heavy'] <= 40) or p['Rings'] < 1 \
        or (no_azole and has_azole(m)) or CATALOG.HasMatch(m)
    if bad: _cache[smi] = None; return None
    sa = sascorer.calculateScore(m)
    if sa > 4.0: _cache[smi] = None; return None
    fp = GEN.GetFingerprint(m)
    mx = float(max(DataStructs.BulkTanimotoSimilarity(fp, TRAIN_FPS)))
    if mx >= 0.50: _cache[smi] = None; return None
    arr = np.zeros(2048, dtype=np.uint8)
    ConvertToNumpyArray(fp, arr)
    d = list(CALC.CalcDescriptors(m))
    d = [MED.iloc[i] if (v != v or abs(v) > 1e10) else v for i, v in enumerate(d)]
    pa = float(MDL.predict_proba(np.hstack([arr, d]).reshape(1, -1))[0, 1])
    scaf = MurckoScaffold.MurckoScaffoldSmiles(mol=m, includeChirality=False)
    r = dict(smiles=smi, prob_active=pa, max_tanimoto=mx, novelty=1 - mx, sa=sa, qed=QED.qed(m),
             scaffold=scaf, scaffold_novel=int(scaf not in TRAIN_SCAF), chemotype=chemotype(smi), **p)
    _cache[smi] = r
    return r

np.random.seed(42)
random.seed(42)

def objvec(r): return np.array([r['prob_active'], r['novelty'], (4.0 - r['sa']) / 4.0, r['qed']])

def nds(F):
    n = len(F); S = [[] for _ in range(n)]; nd = np.zeros(n, int); fronts = [[]]
    for p in range(n):
        gp = np.all(F[p] >= F, axis=1) & np.any(F[p] > F, axis=1)
        gq = np.all(F >= F[p], axis=1) & np.any(F > F[p], axis=1)
        gp[p] = gq[p] = False
        S[p] = list(np.where(gp)[0]); nd[p] = int(gq.sum())
        if nd[p] == 0: fronts[0].append(p)
    i = 0
    while fronts[i]:
        nxt = []
        for p in fronts[i]:
            for q in S[p]:
                nd[q] -= 1
                if nd[q] == 0: nxt.append(q)
        i += 1; fronts.append(nxt)
    return [f for f in fronts if f]

def crowd(F, idx):
    d = np.zeros(len(idx)); A = F[idx]
    for k in range(A.shape[1]):
        o = np.argsort(A[:, k]); d[o[0]] = d[o[-1]] = np.inf
        rg = A[o[-1], k] - A[o[0], k]
        if rg <= 0: continue
        for j in range(1, len(idx) - 1): d[o[j]] += (A[o[j + 1], k] - A[o[j - 1], k]) / rg
    return d

def run_ga(no_azole, side, POP=250, GENS=40, tag='', ckpt=None):
    seeds = [m for m in (Chem.MolFromSmiles(s) for s in act.std_smiles.sample(600, random_state=7)) if m]
    pool = {}; parents = []; t0 = time.time(); tries = 0
    while len(parents) < POP and tries < POP * 150:
        tries += 1
        s = crossover(*random.sample(seeds, 2)) if random.random() < 0.6 else mutate(random.choice(seeds), side)
        if not s: continue
        r = score(s, no_azole)
        if r and s not in pool: pool[s] = r; parents.append(r)
    for g in range(GENS):
        F = np.array([objvec(r) for r in parents]); fr = nds(F)
        rank = np.zeros(len(parents), int)
        for i, f in enumerate(fr):
            for j in f: rank[j] = i
        cd = np.zeros(len(parents))
        for f in fr: cd[f] = crowd(F, f)
        def pick():
            i, j = random.sample(range(len(parents)), 2)
            if rank[i] != rank[j]: return parents[i if rank[i] < rank[j] else j]
            return parents[i if cd[i] > cd[j] else j]
        kids = []; tr_ = 0
        while len(kids) < POP and tr_ < POP * 70:
            tr_ += 1
            if random.random() < 0.65:
                m1 = Chem.MolFromSmiles(pick()['smiles']); m2 = Chem.MolFromSmiles(pick()['smiles'])
                s = crossover(m1, m2) if (m1 and m2) else None
            else:
                m1 = Chem.MolFromSmiles(pick()['smiles']); s = mutate(m1, side) if m1 else None
            if not s or s in pool: continue
            r = score(s, no_azole)
            if r: pool[s] = r; kids.append(r)
        comb = parents + kids; F = np.array([objvec(r) for r in comb]); fr = nds(F); new = []
        for f in fr:
            if len(new) + len(f) <= POP: new += [comb[i] for i in f]
            else:
                o = np.argsort(-crowd(F, f)); new += [comb[f[i]] for i in o[:POP - len(new)]]; break
        parents = new
        if ckpt and (g % 5 == 0 or g == GENS - 1):
            pd.DataFrame(pool.values()).assign(track=tag).to_csv(ckpt, index=False)
    return pool

poolB = run_ga(True, SIDE_NOAZ, tag='azole-free', ckpt='pool_azolefree.csv')
poolA = run_ga(False, SIDE_ALL, tag='unconstrained', ckpt='pool_unconstrained.csv')
dfB = pd.DataFrame(poolB.values()); dfB['track'] = 'Azole-free (primary)'
dfA = pd.DataFrame(poolA.values()); dfA['track'] = 'Unconstrained (control)'
gp = pd.concat([dfB, dfA], ignore_index=True)
gp.to_csv("generated_pool.csv", index=False)

P = np.column_stack([1.0 - gp.prob_active.values, gp.prob_active.values])
grp = np.array([c if c in KEEP else 'POOLED' for c in gp.chemotype])
pv = np.ones((len(gp), 2))
for l in (0, 1):
    a = 1.0 - P[:, l]
    for g in np.unique(grp):
        m = grp == g; ref = CAL.get((g, l))
        if ref is None or len(ref) == 0: ref = CAL[('POOLED', l)]
        pv[m, l] = (len(ref) - np.searchsorted(ref, a[m], side='left') + 1) / (len(ref) + 1)
gp['p_inactive'], gp['p_active'] = pv[:, 0], pv[:, 1]
EPS = 0.05
gp['cp_region'] = np.select(
    [(gp.p_active > EPS) & (gp.p_inactive <= EPS), (gp.p_active <= EPS) & (gp.p_inactive > EPS),
     (gp.p_active > EPS) & (gp.p_inactive > EPS)], ['active', 'inactive', 'both'], default='empty')
gp.to_csv("generated_pool.csv", index=False)

sel = gp[(gp.scaffold_novel == 1) & (gp.max_tanimoto < 0.40) & (gp.sa <= 3.5) & (gp.qed >= 0.50) &
         (gp.prob_active >= 0.80) & (gp.cp_region == 'active')].copy()
sel = sel.sort_values('prob_active', ascending=False).drop_duplicates('scaffold')
B = sel[sel.track == 'Azole-free (primary)'].head(50)
A = sel[sel.track == 'Unconstrained (control)'].head(12)
cand = pd.concat([B, A], ignore_index=True)
cand['cand_id'] = [f"CAL-{i+1:03d}" for i in range(len(cand))]
cand.to_csv("docking_candidates.csv", index=False)

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
    rows.append(dict(smiles=s, max_tanimoto=mx, novelty=1 - mx, sa=sascorer.calculateScore(m),
                     qed=QED.qed(m), scaffold=scaf, scaffold_novel=int(scaf not in TRAIN_SCAF),
                     chemotype=chemotype(s)))
od = pd.DataFrame(rows)
od.to_csv("original_template_enumeration.csv", index=False)

def bits(smis):
    X = []
    for s in smis:
        m = Chem.MolFromSmiles(s)
        a = np.zeros(2048, dtype=np.uint8)
        if m is not None: ConvertToNumpyArray(GEN.GetFingerprint(m), a)
        X.append(a)
    return np.array(X)

tr_act = act.std_smiles.sample(2000, random_state=0).tolist()
B = gp[gp.track == 'Azole-free (primary)'].sample(2000, random_state=0)
A = gp[gp.track == 'Unconstrained (control)'].sample(2000, random_state=0)
groups = [('Training actives', tr_act), ('Original template enumeration', od.smiles.tolist()),
          ('Azole-free track', B.smiles.tolist()), ('Unconstrained track', A.smiles.tolist()),
          ('Selected candidates', cand.smiles.tolist())]
X = np.vstack([bits(g[1]) for g in groups])
lab = np.concatenate([[g[0]] * len(g[1]) for g in groups])
emb = umap.UMAP(n_neighbors=25, min_dist=0.25, metric='jaccard', random_state=42).fit_transform(X)
pd.DataFrame(dict(x=emb[:, 0], y=emb[:, 1], group=lab)).to_csv("chem_space_umap.csv", index=False)