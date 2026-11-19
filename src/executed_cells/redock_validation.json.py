# AUDIT RECORD for `redock_validation.json` - automated reconstruction of the producing cell by the platform lineage tool.
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

import os
import numpy as np
import subprocess
import json
import requests
from rdkit import Chem
from rdkit.Chem import AllChem, rdMolAlign

D = "dock"
os.makedirs(D, exist_ok=True)

# Download 5TZ1 PDB
r = requests.get("https://files.rcsb.org/download/5TZ1.pdb", timeout=60)
with open(os.path.join(D, "5TZ1.pdb"), "wb") as f:
    f.write(r.content)

pdb = open(os.path.join(D, "5TZ1.pdb")).read().splitlines()
prot = [l for l in pdb if l.startswith("ATOM") and l[21] == "A"]
hem = [l for l in pdb if l.startswith("HETATM") and l[17:20] == "HEM" and l[21] == "A"]
vt1 = [l for l in pdb if l.startswith("HETATM") and l[17:20] == "VT1" and l[21] == "A"]

open(os.path.join(D, "receptor_A.pdb"), "w").write("\n".join(prot + hem) + "\nTER\nEND\n")
open(os.path.join(D, "vt1_ref.pdb"), "w").write("\n".join(vt1) + "\nEND\n")

xyz = np.array([[float(l[30:38]), float(l[38:46]), float(l[46:54])] for l in vt1])
ctr = xyz.mean(0)
box = np.maximum(xyz.max(0) - xyz.min(0) + 10.0, 22.0)
np.save(os.path.join(D, "box.npy"), np.vstack([ctr, box]))

subprocess.run(
    ["obabel", os.path.join(D, "receptor_A.pdb"), "-O", os.path.join(D, "receptor.pdbqt"), "-xr", "-p", "7.4"],
    capture_output=True, text=True
)

# Download Vina
vina_url = "https://github.com/ccsb-scripps/AutoDock-Vina/releases/download/v1.2.7/vina_1.2.7_win.exe"
vina_path = os.path.join(D, "vina.exe")
vr = requests.get(vina_url, timeout=120)
with open(vina_path, "wb") as f:
    f.write(vr.content)

# Get VT1 SMILES from RCSB
j = requests.get("https://data.rcsb.org/rest/v1/core/chemcomp/VT1", timeout=60).json()
smi = None
for d in j.get("pdbx_chem_comp_descriptor", []):
    if d["type"] == "SMILES_CANONICAL" and d["program"] == "OpenEye OEToolkits":
        smi = d["descriptor"]
if smi is None:
    for d in j.get("pdbx_chem_comp_descriptor", []):
        if d["type"] == "SMILES":
            smi = d["descriptor"]
            break

open(os.path.join(D, "vt1_smiles.txt"), "w").write(smi)
mh = Chem.AddHs(Chem.MolFromSmiles(smi))
AllChem.EmbedMolecule(mh, randomSeed=42)
AllChem.MMFFOptimizeMolecule(mh)
Chem.MolToMolFile(mh, os.path.join(D, "vt1_gen.sdf"))
subprocess.run(
    ["obabel", os.path.join(D, "vt1_gen.sdf"), "-O", os.path.join(D, "vt1_gen.pdbqt")],
    capture_output=True, text=True
)

import re
import time

rec = os.path.abspath(os.path.join(D, "receptor.pdbqt"))
CTR = ("70.610", "66.284", "4.177")
SZ = ("22.0", "23.5", "22.8")

t = time.time()
result = subprocess.run(
    [os.path.abspath(vina_path), "--receptor", rec,
     "--ligand", os.path.abspath(os.path.join(D, "vt1_gen.pdbqt")),
     "--center_x", CTR[0], "--center_y", CTR[1], "--center_z", CTR[2],
     "--size_x", SZ[0], "--size_y", SZ[1], "--size_z", SZ[2],
     "--exhaustiveness", "32", "--num_modes", "9", "--seed", "42", "--cpu", "4",
     "--out", os.path.abspath(os.path.join(D, "redock_vt1.pdbqt"))],
    capture_output=True, text=True
)

redock_affinity = None
for line in result.stdout.splitlines():
    m = re.match(r"\s*1\s+(-?\d+\.?\d*)", line)
    if m:
        redock_affinity = float(m.group(1))
        break

# Split best mode from vina output
txt = open(os.path.join(D, "redock_vt1.pdbqt")).read().split("MODEL")
best = "MODEL" + txt[1]
open(os.path.join(D, "redock_best.pdbqt"), "w").write(best)
subprocess.run(
    ["obabel", os.path.join(D, "redock_best.pdbqt"), "-O", os.path.join(D, "redock_best.pdb")],
    capture_output=True, text=True
)

tmpl = Chem.MolFromSmiles(smi)
ref = Chem.MolFromPDBFile(os.path.join(D, "vt1_ref.pdb"), removeHs=True, sanitize=False)
pose = Chem.MolFromPDBFile(os.path.join(D, "redock_best.pdb"), removeHs=True, sanitize=False)

def fix(m):
    m = Chem.RemoveHs(m, sanitize=False)
    return AllChem.AssignBondOrdersFromTemplate(tmpl, m)

refb = fix(ref)
poseb = fix(pose)
rmsd = rdMolAlign.CalcRMS(poseb, refb)

hem_lines = [l for l in open(os.path.join(D, "5TZ1.pdb")) if l.startswith("HETATM") and l[17:20] == "HEM" and l[21] == "A"]
fe_line = [l for l in hem_lines if l[12:16].strip() == "FE"][0]
fexyz = np.array([float(fe_line[30:38]), float(fe_line[38:46]), float(fe_line[46:54])])

def tet_n_to_fe(m, fe):
    conf = m.GetConformer()
    patt = Chem.MolFromSmarts("n1cnnn1")
    mt = m.GetSubstructMatches(patt)
    best_d = 1e9
    for mm in mt:
        for a in mm:
            p = conf.GetAtomPosition(a)
            d = np.linalg.norm(np.array([p.x, p.y, p.z]) - fe)
            best_d = min(best_d, d)
    return best_d

crystal_N_Fe = tet_n_to_fe(refb, fexyz)
docked_N_Fe = tet_n_to_fe(poseb, fexyz)

json.dump(
    dict(
        redock_rmsd=round(float(rmsd), 3),
        redock_affinity=redock_affinity,
        crystal_N_Fe=round(crystal_N_Fe, 2),
        docked_N_Fe=round(docked_N_Fe, 2),
        pdb='5TZ1',
        resolution=2.0
    ),
    open("redock_validation.json", "w"),
    indent=2
)