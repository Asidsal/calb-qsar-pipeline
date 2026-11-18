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

## Research files

| Path | Contents |
| --- | --- |
| `data/` | Curated MIC data, training/test partitions and chemotype information |
| `notebooks/` | Data curation, feature construction and model training |
| `models/` | Serialized predictive assets where included |
| `src/` | Molecular generation and documented analysis code |
| `results/` | Validation metrics, molecular-generation and candidate data |
| `docking/` | Protein inputs, validated docking box and reference structures |
| `md/` | Prepared CYP51–ligand structures, 50 ns analyses and MM/GBSA results |
| `docs/retrosynthesis/` | CAL-028 synthetic planning: scheme, building blocks and route steps |

Detailed MM/GBSA summaries are in [`md/mmgbsa/results/`](md/mmgbsa/results/). The complete production DCD trajectories (8.44 GiB across CAL-028 and CAL-015) and full per-frame energy logs are kept in separate project storage; lighter-derived datasets, visualization files and reproducible workflows can be kept in GitHub.

## Reproduction

Use the versions recorded by the workflow; for the core QSAR pipeline, run the supplied notebooks in order after installing the indicated dependencies. The GA and result files use the recorded random seed and structural-filter definition. The MD analyses read AMBER `.prmtop` and DCD inputs. MM/GBSA uses the AmberTools `MMPBSA.py` program, with the heme included as receptor cofactor.

The proposed retrosynthesis is supplied as a chemical-design route and intermediate structure set. The route is separate from the calculated binding-mode and QSAR results.

## Citation and license

This is an evolving research code/data repository; a publication venue and acceptance status are **not specified**. Code is released under the MIT license where covered by [LICENSE](LICENSE). ChEMBL-derived data retain their applicable terms and attribution.
