# CYP51 single-trajectory MM/GBSA

The CAL-028 and CAL-015 systems were evaluated independently using AmberTools 26 `MMPBSA.py`, sampling 200 frames throughout the corresponding 50 ns production trajectory (five 10 ns blocks, 40 frames each).

## Configuration

- GB-OBC2 (`igb=5`), `mbondi2` intrinsic radii, `saltcon=0.150` M.
- Single-trajectory complex/receptor/ligand protocol, stripping explicit solvent and ions from the calculation snapshots.
- Receptor includes CYP51, its heme cofactor and the bonded proximal cysteine; CAL-028 is retained as its modeled anionic carboxylate.
- Components: `VDWAALS`, `EEL`, `EGB` and `ESURF`. `DELTA TOTAL` is the reported MM/GBSA endpoint (conformational entropy not included).

## Completed 50 ns results

| Compound | N frames | MM/GBSA ΔTOTAL, kcal mol⁻¹ |
| --- | ---: | ---: |
| CAL-028 | 200 | −60.42092 ± 3.72364 |
| CAL-015 | 200 | −59.46166 ± 3.44685 |

[Full energy component summary](results/MMGBSA_both_compounds_summary.csv) · [10 ns block summary](results/MMGBSA_10ns_block_summary.csv).

The energy profiles complement the hydrogen-bond and protein-contact analyses. Favorable van der Waals and nonpolar terms account for much of the net energy for both ligands; polar solvation offsets a substantial part of the electrostatic contribution.

The source DCD trajectories, AMBER parameter/topology records, per-frame AmberTools outputs and complete calculation logs are retained in the project's separate storage. The tables above report the executed 50 ns results.

## Computational reference

Miller, B. R. III et al. *MMPBSA.py: An Efficient Program for End-State Free Energy Calculations*. J. Chem. Theory Comput. **2012**, 8, 3314–3321. DOI: 10.1021/ct300418h.
