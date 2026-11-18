# Molecular dynamics analysis — CYP51 with CAL-028 and CAL-015

Independent 50 ns production trajectories were analyzed for each CYP51–ligand complex.

## Summary

| Measurement | CAL-028 | CAL-015 |
| --- | ---: | ---: |
| Ligand RMSD, first production frame (Å) | 0.919 ± 0.247 | 1.579 ± 0.185 |
| Ligand RMSD, original docked pose (Å) | 2.158 ± 0.155 | 2.524 ± 0.164 |
| Protein Cα RMSD (Å) | 1.389 ± 0.185 | 1.513 ± 0.142 |
| Protein Cα radius of gyration (Å) | 22.571 ± 0.067 | 22.581 ± 0.065 |
| Bound-ligand SASA (Å²) | 32.91 ± 7.03 | 47.47 ± 14.05 |
| Direct protein–ligand H bonds per sampled frame | 3.538 | 1.115 |
| Proximal heme Fe–S distance (Å) | 2.644 ± 0.061 | 2.658 ± 0.067 |

The full numerical dataset is [`derived_summary.csv`](derived_summary.csv). Both protein trajectories contain 484 analyzed Cα positions. Time-series analysis used 5,000 frames per compound; solvent exposure was sampled from 200 frames and geometric H-bond analysis from 1,000 frames.

## Individual binding modes

**CAL-028.** Its carboxylate-centered interaction network involves Tyr118, Ser378, Ile379 and Phe380. The ligand retains a localized pocket pose while changing its position relative to heme iron.

**CAL-015.** A persistent Met508 backbone hydrogen bond accompanies extensive contacts with residues lining the CYP51 cavity. The stereodefined simulated structure retains its S stereocenter and E-fluoroalkene geometry.

The original AMBER topology and full DCD source files are retained in separate project storage; the lightweight MD metric tables are released here. See [MM/GBSA methods and results](../mmgbsa/README.md).
