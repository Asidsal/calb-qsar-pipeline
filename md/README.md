# MD-ready complexes

Docked protein-ligand starting structures for molecular dynamics.

- `complexes/<id>_complex.pdb` - chain A of PDB 5TZ1 (C. albicans CYP51) + heme (HEM) + docked ligand (resname LIG,
  chain L, residue 901), best AutoDock Vina pose, heavy atoms only. Real element symbols; CONECT records for the ligand
  and the heme.
- `ligands/<id>_ligand_pose.sdf` - the same pose as a stand-alone ligand with bond orders and explicit hydrogens
  (neutral form; assign the ionisation state for the simulated pH yourself).
- `docked_smiles.csv` - SMILES of each docked pose (stereo as docked) next to the GA-generated SMILES.
- `MANIFEST.sha256` - SHA-256 of every file in the original package.

Ids CAL-028/015/032/017 are Tier 1; CAL-053/058/051 Tier 2; CAL-025/050 carry structural alerts; fluconazole and
posaconazole are reference azoles. Original DCD trajectories and AMBER topologies are stored externally; derived 50 ns results, analysis scripts and visualization outputs are now included under `md/analysis_50ns/`, `md/workflows/`, `md/visualization_3d/` and `md/mmgbsa/`.

## Stereochemical identities of docked poses
No specified stereocentre was inverted by docking. But the GA wrote SMILES that leave some stereo
undefined, and docking evaluated one arbitrary geometry for each:

- **CAL-015 and CAL-032**: the fluoroalkene `CCC=C(F)` has no E/Z specified in the generated SMILES.
  Docking embedded one geometry for each and they are *opposite* ones: CAL-015 was docked as
  `CC/C=C(/F)[C@H](CNC(=O)c1ccn(CC2CCNCC2)n1)NC(=O)c1ccc(Cl)cc1` and CAL-032 as `CC/C=C(\F)[C@H](CNC(=O)c1ccn(CC2CCNCC2)n1)NC1CCCCCCC1`. The two poses are therefore not
  like-for-like.
- **CAL-017**: the alkene `C=C(C)` and the piperidin-3-yl carbon were both unspecified. Docked as
  `CN[C@H](C)CNC(=O)c1ccn(Cc2ccc(/C=C(\C)[C@@H]3CCCNC3)cc2)n1`.
- **CAL-028** has one defined stereocentre and no ambiguous alkene.

For CAL-015, CAL-032 and CAL-017, simulate **both** alkene isomers (and both piperidine epimers for
CAL-017) before drawing conclusions, and specify the isomer in any retrosynthesis.

## Heme geometry of the docked poses (angstrom)
Genuine azole coordination has a ligand N within about 2.1-2.3 A of the iron.

| compound | closest atom to Fe | nearest pyridine-like aromatic N | nearest any N |
|---|---|---|---|
| CAL-028 | F1 at 3.29 | 7.91 | 7.49 |
| CAL-015 | F1 at 3.20 | 9.32 | 5.12 |
| CAL-032 | C15 at 3.05 | 12.39 | 7.74 |
| CAL-017 | N1 at 6.34 | 13.07 | 6.34 |
| CAL-025 | C19 at 3.22 | 4.04 | 4.04 |
| fluconazole | F2 at 4.38 | 5.79 | 5.39 |
| posaconazole | C20 at 4.74 | 5.29 | 5.29 |

**None of the four Tier 1 poses puts a nitrogen anywhere near the iron** (5.1 A at best), and the
closest atoms are fluorine or carbon at about 3.0-3.3 A, which is a contact and not coordination. Vina
also fails to coordinate fluconazole (nearest triazole N 5.79 A) and posaconazole (5.29 A), although it
does reproduce the VT-1161 tetrazole (2.60 A) in the redocking control. Vina has no term for the Fe-N
coordinate bond, so these distances say nothing for or against a heme-binding mechanism. They do mean
the imidazole/triazole-excluded candidates cannot be assumed to inhibit CYP51 by the azole route. Completed 50 ns MD for CAL-028 and the simulated E isomer of CAL-015 is documented under `analysis_50ns/`; the other candidates have no corresponding MD analysis in this repository.
