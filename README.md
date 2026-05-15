# Scaffold-Aware QSAR Pipeline for *Candida albicans* Antifungal Activity

This repository contains the complete computational pipeline accompanying the manuscript:

> **Scaffold-Aware QSAR Modeling of *Candida albicans* Antifungal Activity: Applicability Domain Characterization, Model Validation, and Genetic Algorithm-Guided Design of Novel Triazole Candidates**  
> Abubakar Siddiq Salihu and Wan Mohd Nuzul Hakimi Wan Salleh  
> *Journal of Chemical Information and Modeling*, 2026 (under review)

The pipeline curates 11,245 MIC records from ChEMBL 36, trains and benchmarks seven machine learning classifiers under a stratified Murcko scaffold split, characterises the model applicability domain (AD) using a Tanimoto k-nearest-neighbour index, and applies a pharmacophore-constrained genetic algorithm to generate 15 novel triazole antifungal candidates all within the model AD (predicted probability 0.709-0.959).

---

## Repository structure

```
calb-qsar-pipeline/
├── notebooks/
│   ├── 01_data_curation.ipynb           ChEMBL 36 extraction and curation
│   ├── 02_scaffold_split_features.ipynb Stratified scaffold split + ECFP4/descriptor generation
│   ├── 03_model_training.ipynb          7-model benchmark + 5-fold CV + scaffold test evaluation
│   ├── 04_shap_ad_analysis.ipynb        SHAP feature importance + AD characterisation
│   └── 05_GA_design_colab.ipynb         GA-guided de novo triazole generation (Google Colab)
├── data/
│   └── CAlbicans_curated.csv            Curated dataset (11,245 compounds, binary MIC labels)
├── requirements.txt
└── README.md
```

---

## Requirements

```bash
pip install -r requirements.txt
```

Python 3.10+ recommended. Notebook 5 is designed for Google Colab (CPU runtime).

---

## Data

ChEMBL 36 is required to reproduce the curation step (Notebook 1). Download the SQLite version from the [ChEMBL downloads page](https://chembl.gitbook.io/chembl-interface-documentation/downloads) and update `DB_PATH` in Notebook 1 accordingly. The curated output `data/CAlbicans_curated.csv` is provided directly so Notebooks 2-4 can be run without a local ChEMBL installation.

---

## Usage

Run notebooks in order: 01 -> 02 -> 03 -> 04. Notebook 05 is standalone and runs on Google Colab; upload the four files listed in its header cell when prompted.

---

## Contact

Abubakar Siddiq Salihu - abusiddiq.salihu@umyu.edu.ng  
ORCID: 0000-0002-4425-7524  
Scopus ID: 58105245900
