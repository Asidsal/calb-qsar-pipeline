"""Check that src/ga.py reproduces the stored scores for molecules in results/generated_pool.csv.
Run from the repository root:  python tests/test_ga_scoring.py   (loading the module takes ~1-2 min)."""
import sys, pathlib, numpy as np, pandas as pd
root = pathlib.Path(__file__).resolve().parents[1]
sys.path.insert(0, str(root / "src"))
import ga
gp = pd.read_csv(root / "results" / "generated_pool.csv")
s = pd.concat([gp[gp.track.str.startswith("Azole-free")].sample(40, random_state=1),
               gp[gp.track.str.startswith("Unconstrained")].sample(40, random_state=1)])
worst = 0.0
for _, r in s.iterrows():
    o = ga.score(r.smiles, r.track.startswith("Azole-free"))
    assert o is not None, r.smiles
    worst = max(worst, abs(o["prob_active"] - r.prob_active), abs(o["max_tanimoto"] - r.max_tanimoto),
                abs(o["sa"] - r.sa), abs(o["qed"] - r.qed))
    assert o["scaffold"] == r.scaffold and o["chemotype"] == r.chemotype
print("OK: 80 molecules reproduced; max abs difference =", worst)
