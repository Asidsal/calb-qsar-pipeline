# Result-to-code index

This index maps reported data or result files to the Python analysis code retained in `src/executed_cells/` and lists their recorded inputs. It is intended to help locate and execute the relevant computational steps. The core curation, feature preparation and model-training workflow is provided separately in `notebooks/01_data_curation.ipynb`, `notebooks/02_scaffold_split_features.ipynb`, and `notebooks/03_model_training.ipynb`.

| Result file | Inputs recorded | Code file |
|---|---|---|
| `chemotype_assignments.csv` | CAlbicans_curated.csv | `chemotype_assignments.csv.py` |
| `generalization_benchmark.csv` | chemotype_assignments.csv, generalization_sens.csv, _gb_partial.csv | `generalization_benchmark.csv.py` |
| `cluster_split_repeats.csv` | pool_meta.csv, y_all.npy | `cluster_split_repeats.csv.py` |
| `conformal_results.csv` | chemotype_assignments.csv, X_train_ecfp4.npy, X_test_ecfp4.npy, X_train_desc.csv, X_test_desc.csv, y_train.npy, y_test.npy, train.csv, test.csv | `conformal_results.csv.py` |
| `conformal_chemotype_validity.csv` | chemotype_assignments.csv, cal_idx.npy, prop_idx.npy | `conformal_chemotype_validity.csv.py` |
| `conformal_vs_knnAD.csv` | chemotype_assignments.csv, cal_idx.npy, prop_idx.npy | `conformal_vs_knnAD.csv.py` |
| `conformal_predictor.pkl` | chemotype_assignments.csv | `conformal_predictor.pkl.py` |
| `chemotype_only_baseline.csv` | chemotype_assignments.csv | `chemotype_only_baseline.csv.py` |
| `external_panel_v2.csv` | conformal_predictor.pkl, partA_antifungal_validation.csv, train.csv, test.csv | `external_panel_v2.csv.py` |
| `external_validation_v2.csv` | conformal_predictor.pkl, chemotype_assignments.csv, X_train_desc.csv | `external_validation_v2.csv.py` |
| `erg11_curated.csv` | chembl_erg11_raw.pkl, 01_data_curation.ipynb | `erg11_curated.csv.py` |
| `human_cyp51_curated.csv` | 01_data_curation.ipynb, chembl_hcyp51_raw.pkl | `human_cyp51_curated.csv.py` |
| `docking_scores_enrichment.csv` | chemotype_assignments.csv, conformal_predictor.pkl | `docking_scores_enrichment.csv.py` |
| `docking_validation.csv` | chemotype_assignments.csv, conformal_predictor.pkl | `docking_validation.csv.py` |
| `redock_validation.json` | chemotype_assignments.csv, conformal_predictor.pkl | `redock_validation.json.py` |
| `candidates_with_docking.csv` | - | `candidates_with_docking.csv.py` |
| `docking_candidates.csv` | chemotype_assignments.csv, conformal_predictor.pkl | `docking_candidates.csv.py` |
| `candidates_scored.csv` | chemotype_assignments.csv, conformal_predictor.pkl | `candidates_scored.csv.py` |
| `final_candidates.csv` | chemotype_assignments.csv, conformal_predictor.pkl, external_panel_v2.csv, fig01_baseline_audit.png, fig02_dual_axis_dataset.png, fig03_conformal_validity.png, fig04_generalization.png, fig05_external_validation.png, fig06_docking_enrichment.png, fig07_ga_convergence.png, fig08_chemical_space.png, fig09_candidate_structures.png, fig10_docking_poses.png, generalization_benchmark.csv, conformal_vs_knnAD.csv, conformal_results.csv, conformal_chemotype_validity.csv, external_validation_v2.csv, docking_validation.csv, README.md, audit_report.md, dual_axis_report.md, conformal_report.md, generalization_report.md, external_validation_report.md, docking_report.md, candidate_dossiers.md, md_inputs_README.md, final_candidates.csv, candidates_scored.csv, cluster_split_repeats.csv, baseline_metrics.csv, external_validation_summary.csv, docking_scores_enrichment.csv, redock_validation.json, pose_contacts.csv, original_template_enumeration.csv, erg11_curated.csv, human_cyp51_curated.csv, ga_history.csv, manuscript_outline.md, internal_notes.md, figure_captions.md, manuscript_tables.xlsx, generated_pool.csv, chem_space_umap.csv, docking_candidates.csv, pool_accounting.json, RESULTS_LEDGER.csv, draft.txt, manuscript_tables.xlsx, F7_chemical_space.png, F6_ga_convergence.png, F4_external_validation.png, README_repository.md, README.md, generalization_sens.csv, _gb_partial.csv, chembl_hcyp51_raw.pkl, chembl_erg11_raw.pkl, test_knn_sim.npy, repro_svm.json, repro_test_metrics.json | `final_candidates.csv.py` |
| `dual_axis_overlap.csv` | chembl_erg11_raw.pkl, CAlbicans_curated.csv | `dual_axis_overlap.csv.py` |
| `pose_contacts.csv` | - | `pose_contacts.csv.py` |
| `chem_space_umap.csv` | chemotype_assignments.csv, conformal_predictor.pkl | `chem_space_umap.csv.py` |
| `original_template_enumeration.csv` | chemotype_assignments.csv, conformal_predictor.pkl | `original_template_enumeration.csv.py` |
| `baseline_metrics.csv` | repro_svm.json, repro_test_metrics.json, model_results.json | `baseline_metrics.csv.py` |
| `pool_accounting.json` | chemotype_assignments.csv, Whole-cell antifungal QSAR models learn chemotype, not activity scope characteri.docx | `pool_accounting.json.py` |

Execution environment, datasets and validation commands are described in the project README. Outputs are preserved in `results/` and `data/`.
