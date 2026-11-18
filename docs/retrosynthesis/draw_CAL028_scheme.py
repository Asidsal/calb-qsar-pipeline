#!/usr/bin/env python3
"""Draw a chemical-structure scheme for the proposed convergent CAL-028 route.

Dependencies: rdkit, matplotlib, pillow.
Run: python docs/retrosynthesis/draw_CAL028_scheme.py
Output: docs/retrosynthesis/CAL028_proposed_scheme.png
"""
import csv
from pathlib import Path
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.offsetbox import OffsetImage, AnnotationBbox
from rdkit import Chem
from rdkit.Chem import Draw

ROOT = Path(__file__).resolve().parent
with (ROOT / "CAL028_intermediates.csv").open(newline="", encoding="utf-8") as f:
    MOL = {row["id"]: Chem.MolFromSmiles(row["isomeric_smiles"])
           for row in csv.DictReader(f)}
for code, mol in MOL.items():
    if mol is None:
        raise ValueError(f"Invalid SMILES in {code}")

fig, ax = plt.subplots(figsize=(16, 11.8))
ax.set_xlim(0, 1)
ax.set_ylim(0, 1)
ax.axis("off")
fig.patch.set_facecolor("white")

def structure(code, x, y, scale, label):
    img = Draw.MolToImage(MOL[code], size=(380, 260))
    ax.add_artist(AnnotationBbox(OffsetImage(img, zoom=scale), (x, y),
                  frameon=False, xycoords="axes fraction"))
    ax.text(x, y - .100, label, ha="center", va="top", size=9.5,
            weight="bold", transform=ax.transAxes)

def reaction(x1, x2, y, reagent, fontsize=8):
    ax.annotate("", xy=(x2, y), xytext=(x1, y), xycoords="axes fraction",
                arrowprops={"arrowstyle": "-|>", "lw": 1.2, "color": "black"})
    ax.text((x1+x2)/2, y+.023, reagent, ha="center", va="bottom",
            size=fontsize, transform=ax.transAxes)

ax.text(.02,.985,"PROPOSED CONVERGENT SYNTHESIS OF CAL-028",size=17,
        weight="bold",va="top",transform=ax.transAxes)
ax.text(.02,.945,"(A) N-benzylpyrazole acid branch",size=12,weight="bold",
        va="top",transform=ax.transAxes)
for item in [("A0",.14,"A0 · pyrazole ester"),
             ("A1",.50,"A1 · N-benzyl ester"),
             ("A2",.86,"A2 · pyrazole acid")]:
    structure(item[0], item[1], .82, .35, item[2])
reaction(.27,.37,.825,"Difluorobenzyl bromide / K₂CO₃")
reaction(.63,.74,.825,"LiOH / THF–H₂O")
ax.plot([.02,.98],[.68,.68],color=".8",lw=.7,transform=ax.transAxes)
ax.text(.02,.665,"(B) Stereodefined Nα-methyl-L-Dap branch",size=12,
        weight="bold",va="top",transform=ax.transAxes)
for code,x,label in [("B0",.105,"B0 · Fmoc/Boc precursor"),
                     ("B1",.302,"B1 · methyl ester"),
                     ("B2",.499,"B2 · Nα-Me amine"),
                     ("B3",.696,"B3 · pyrazine amide"),
                     ("B4",.895,"B4 · Nβ amine")]:
    structure(code,x,.535,.245,label)
for x1,x2,t in [(.185,.22,"Esterify"),(.381,.418,"Fmoc off"),
                (.58,.618,"Pyrazine acid / coupling"),(.774,.81,"Boc off")]:
    reaction(x1,x2,.535,t,6.5)
ax.plot([.02,.98],[.39,.39],color=".8",lw=.7,transform=ax.transAxes)
ax.text(.02,.372,"(C) Convergent amide assembly and ester hydrolysis",size=12,
        weight="bold",va="top",transform=ax.transAxes)
structure("A2",.13,.203,.29,"A2 · acid")
ax.text(.267,.203,"+",size=20,ha="center",va="center",transform=ax.transAxes)
structure("B4",.36,.203,.29,"B4 · amine")
reaction(.47,.54,.204,"Amide coupling")
structure("B5",.67,.203,.32,"B5 · methyl ester")
reaction(.785,.82,.204,"LiOH")
structure("B6",.91,.203,.31,"(S)-CAL-028")
ax.text(.02,.016,"Fmoc = fluorenylmethoxycarbonyl; Boc = tert-butoxycarbonyl; "
        "Dap = diaminopropionic acid. Proposed route.",size=8,
        transform=ax.transAxes)
fig.savefig(ROOT / "CAL028_proposed_scheme.png",dpi=250,
            bbox_inches="tight",facecolor="white")
print(ROOT / "CAL028_proposed_scheme.png")
