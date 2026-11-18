# Chemotype-aware antifungal QSAR, molecular design and CYP51 binding studies

A reproducible computational discovery workflow for *Candida albicans* whole-cell antifungal activity: curated ChEMBL activity records, scaffold- and chemotype-aware QSAR validation, conformal applicability assessment, molecular generation, structure-based prioritization, and 50 ns CYP51 molecular dynamics and MM/GBSA characterization of CAL-028 and CAL-015.

This repository contains research **code, structured input data, model/result tables, molecular structures and analysis workflows**. The manuscript and its Supporting Information are maintained separately and are **not included here**.

## Research scope and findings

- **Data:** 11,245 curated ChEMBL 36 MIC records (4,416 active; 6,829 inactive), activity cutoff 8 μg/mL.
- **Generalization:** scaffold-split MCC 0.667; leave-one-chemotype-out tests examine the effect of chemical distribution shift and define the model's transferability across chemotypes.
- **Uncertainty assessment:** Mondrian conformal prediction and nearest-neighbor applicability assessment, including subgroup-conditional coverage analyses.
- **Design:** graph genetic algorithm with BRICS-fragment mutation, multiobjective survivor selection and two molecular-generation tracks. The primary structural filter excludes imidazole and triazole motifs; other heterocycles, including pyrazoles, are permitted.
- **Docking:** CYP51 structure 5TZ1, reference redocking and enrichment validation, followed by candidate prioritization.
- **Dynamics:** CAL-028 and CAL-015 were individually simulated for **50 ns**, with RMSD, RMSF, SASA, hydrogen-bond timelines, residue-contact maps, heme geometry and pose-evolution analyses.
- **MM/GBSA:** 200 frames per compound over five 10 ns blocks, GB-OBC2 (igb=5), mbondi2 radii and 0.150 M salt. Frame-averaged ΔTOTAL: **CAL-028 −60.42 ± 3.72 kcal/mol** and **CAL-015 −59.46 ± 3.45 kcal/mol**, reported as MM/GBSA energies without entropic correction.
- **Synthesis planning:** a convergent **proposed retrosynthetic route for CAL-028**, with a stereodefined N-methyl-L-Dap building block, N-benzylpyrazole coupling fragment and amide-convergence strategy. See [proposed synthesis scheme](docs/retrosynthesis/CAL028_proposed_scheme.md).

## Current public repository layout

| Path | Public contents |
| --- | --- |
| `data/CAlbicans_curated.csv` | Curated ChEMBL whole-cell MIC records |
| `notebooks/` | Curation, split/feature development, classification and the earlier GA design notebook |
| `docs/retrosynthesis/` | Proposed CAL-028 convergent synthesis scheme, stereodefined intermediate SMILES, route-step CSV and RDKit drawing script |
| `md/analysis_50ns/` | Summary of the independent CAL-028 and CAL-015 50 ns trajectory measurements |
| `md/mmgbsa/` | Completed MM/GBSA component means and five time-block results |
| `requirements.txt` | Core notebook environment specification |

The extended research archive also preserves the molecular-generation source, detailed docking records, trained model, complete derived trajectories and visualization resources. Those larger research assets are distributed separately from the current public GitHub files. The full production DCD/AMBER files and per-frame GB results remain in separate project storage.

## Research reproducibility and interpretation

The repository supplies primary curation and modelling notebooks and publishes the exact derived MD/MM/GBSA statistics above. The simulation sources use AMBER topology and DCD trajectories; MM/GBSA uses the AmberTools `MMPBSA.py` single-trajectory protocol with heme included in the receptor. The proposed CAL-028 retrosynthesis records specific intermediate structures and transformations, with a structure-bearing diagram reproducible from the RDKit drawing script.

The computational findings are interpreted in the QSAR, docking, MD and MM/GBSA frameworks described in their source records; their role is to prioritize and characterize new antifungal chemical matter.

## Citation and license

This is an evolving research code/data repository; a publication venue and acceptance status are **not specified**. Code is released under the MIT license where covered by [LICENSE](LICENSE). ChEMBL-derived data retain their applicable terms and attribution.

## Source evidence for synthesis planning and molecular dynamics

- [ASKCOS original CAL-028 retrosynthetic tree export](docs/retrosynthesis/ASKCOS_CAL028_treeResults.json): 200 trees from the [shared ASKCOS result](https://askcos.mit.edu/results?shared=0a31ab2d-cd9b-4b72-a129-7f43a8fd4bfb), with [reproducible route-count audit](docs/retrosynthesis/ASKCOS_CAL028_tree_audit.md) (56 four-step, 144 five-step, 22 distinct reaction transformations; 68/66/66 first-disconnection co-reactants).
- [Executed 50 ns OpenMM production protocol](md/protocol/README.md) and its original [CAL-028](md/protocol/CAL028_50ns.py) and [CAL-015](md/protocol/CAL015_production_50ns_1fs.py) scripts. These document 2 fs and 1 fs production steps, respectively, 300 K and 1 bar Langevin/Monte Carlo conditions, PME and 1 nm cutoff, explicit TIP3P water and the periodic truncated-octahedral simulation cell.

The manuscript and Supporting Information are **not** part of this repository.
