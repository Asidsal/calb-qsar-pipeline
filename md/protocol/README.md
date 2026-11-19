<h1 align="center">50 ns CYP51 molecular dynamics protocol</h1>

This directory contains the OpenMM production scripts, simulation settings and input-state records for the CAL-028 and CAL-015 CYP51 complexes.


## Contents

- [Production settings](#production-settings)
- [Production scripts](#production-scripts)
- [Protein force field](#protein-force-field)
- [Simulation data](#simulation-data)

## Production settings

| Parameter | CAL-028 | CAL-015 |
| --- | --- | --- |
| Simulation engine | OpenMM (CAL-028 manifest: 8.6.1) | OpenMM |
| Protein force field | AMBER ff14SB | AMBER ff14SB |
| Production trajectory | 50 ns | 50 ns |
| Timestep | 2 fs | 1 fs |
| Steps per 1 ns segment | 500,000 | 1,000,000 |
| Ensemble | NPT | NPT |
| Integrator | LangevinMiddleIntegrator | LangevinMiddleIntegrator |
| Temperature | 300 K | 300 K |
| Langevin friction | 1 ps⁻¹ | 1 ps⁻¹ |
| Pressure | 1 bar | 1 bar |
| Pressure coupling | MonteCarloBarostat, interval 25 steps | MonteCarloBarostat, interval 25 steps |
| Electrostatics | PME | PME |
| Nonbonded cutoff | 1.0 nm | 1.0 nm |
| Constraints | Bonds involving hydrogen (`HBonds`) | `HBonds` |
| GPU | CUDA mixed precision | CUDA mixed precision |
| DCD output | Every 5,000 steps (10 ps) | Every 10,000 steps (10 ps) |
| Samples over 50 ns | 5,000 | 5,000 |
| Atom count | 75,277 | 75,288 |
| Water molecules | 22,451 | 22,451 |
| Nominal starting box | 102.730774 Å sides, 109.471219° angles | Same |
| Periodic cell | Truncated octahedron | Truncated octahedron |
| Water | Explicit three-site TIP3P | Explicit three-site TIP3P |
| Ligand net charge | −1 e | +1 e |
| Heme | IC6 ferric, five-coordinate bonded Fe–S and Fe–N | IC6 ferric, same connectivity |

The starting box is extracted from the saved AMBER topologies; dimensions fluctuate under NPT dynamics. Water oxygen charge −0.834 e, each water hydrogen +0.417 e. Counterions neutralize the solvated complexes.

## Production scripts

- [`CAL028_50ns.py`](CAL028_50ns.py): CAL-028 production script, 2 fs, segmented GPU execution, topology and input-state SHA-256 validation
- [`CAL015_production_50ns_1fs.py`](CAL015_production_50ns_1fs.py): CAL-015 production script, 1 fs timestep, segmented checkpoints and heme-geometry monitoring
- [`CAL028_protocol.json`](CAL028_protocol.json): simulation manifest, including OpenMM version, timestep, and hashes

CAL-028 used a 10 ps post-equilibration warm-up excluded from the production timeline. CAL-015 production started from the endpoint of 1 fs NPT equilibration. Both systems used AMBER `.prmtop` topology files containing the force-field parameters.

## Protein force field

The CYP51 protein was parameterized using the AMBER ff14SB force field (Maier et al., 2015; DOI: [10.1021/acs.jctc.5b00255](https://doi.org/10.1021/acs.jctc.5b00255)). Ligands were modeled with GAFF2 and AM1-BCC charges, and the ferric heme group with IC6 parameters. Both production AMBER topologies contain 484 protein Cα atoms assigned atom type `CX`, with no CMAP correction entries.

## Simulation data

The ~8.44 GiB DCD production trajectories and full AMBER parameter/coordinate and checkpoint files are retained in project storage. This repository contains the processing outputs and analysis scripts.
