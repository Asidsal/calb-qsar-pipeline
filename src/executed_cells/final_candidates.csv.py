# AUDIT RECORD for `final_candidates.csv` - automated reconstruction of the producing cell by the platform lineage tool.
# Absolute local paths were replaced by P(<filename>); nothing else was edited. This file has NOT been re-run.
# Inputs recorded: chemotype_assignments.csv, conformal_predictor.pkl, external_panel_v2.csv, fig01_baseline_audit.png, fig02_dual_axis_dataset.png, fig03_conformal_validity.png, fig04_generalization.png, fig05_external_validation.png, fig06_docking_enrichment.png, fig07_ga_convergence.png, fig08_chemical_space.png, fig09_candidate_structures.png, fig10_docking_poses.png, generalization_benchmark.csv, conformal_vs_knnAD.csv, conformal_results.csv, conformal_chemotype_validity.csv, external_validation_v2.csv, docking_validation.csv, README.md, audit_report.md, dual_axis_report.md, conformal_report.md, generalization_report.md, external_validation_report.md, docking_report.md, candidate_dossiers.md, md_inputs_README.md, final_candidates.csv, candidates_scored.csv, cluster_split_repeats.csv, baseline_metrics.csv, external_validation_summary.csv, docking_scores_enrichment.csv, redock_validation.json, pose_contacts.csv, original_template_enumeration.csv, erg11_curated.csv, human_cyp51_curated.csv, ga_history.csv, manuscript_outline.md, internal_notes.md, figure_captions.md, manuscript_tables.xlsx, generated_pool.csv, chem_space_umap.csv, docking_candidates.csv, pool_accounting.json, RESULTS_LEDGER.csv, draft.txt, manuscript_tables.xlsx, F7_chemical_space.png, F6_ga_convergence.png, F4_external_validation.png, README_repository.md, README.md, generalization_sens.csv, _gb_partial.csv, chembl_hcyp51_raw.pkl, chembl_erg11_raw.pkl, test_knn_sim.npy, repro_svm.json, repro_test_metrics.json
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

import pandas as pd
import numpy as np
import pickle
from rdkit import Chem, DataStructs
from rdkit.Chem import Descriptors, rdMolDescriptors, Crippen, AllChem
from rdkit.Chem.FilterCatalog import FilterCatalog, FilterCatalogParams

# Load inputs
m = pd.read_csv(P("candidates_scored.csv"))
panel = pd.read_csv(P("external_panel_v2.csv"))

fin = m[m.docking_supported == 1].sort_values(['selectivity', 'affinity']).reset_index(drop=True)

p2 = FilterCatalogParams()
for c in (FilterCatalogParams.FilterCatalogs.NIH, FilterCatalogParams.FilterCatalogs.CHEMBL):
    p2.AddCatalog(c)
CAT2 = FilterCatalog(p2)

GEN = AllChem.GetMorganGenerator(radius=2, fpSize=2048)

af = panel[panel.panel == 'Antifungal']
af_fps = [(r['query'], GEN.GetFingerprint(Chem.MolFromSmiles(r.std_smiles)))
          for _, r in af.iterrows()
          if Chem.MolFromSmiles(r.std_smiles)]

rows = []
for _, r in fin.iterrows():
    mol = Chem.MolFromSmiles(r.smiles)
    fp = GEN.GetFingerprint(mol)
    nn, ns = max(((n, DataStructs.TanimotoSimilarity(fp, f)) for n, f in af_fps), key=lambda x: x[1])
    lip = sum([Descriptors.MolWt(mol) > 500, Crippen.MolLogP(mol) > 5,
               rdMolDescriptors.CalcNumHBD(mol) > 5, rdMolDescriptors.CalcNumHBA(mol) > 10])
    alerts = [e.GetDescription() for e in CAT2.GetMatches(mol)]
    rows.append(dict(
        cand_id=r.cand_id, track=r.track, smiles=r.smiles,
        prob_active=r.prob_active, p_active=r.p_active,
        affinity_Calbicans=r.affinity, affinity_human=r.affinity_human, selectivity=r.selectivity,
        max_tanimoto_training=r.max_tanimoto, nearest_antifungal=nn, nearest_antifungal_Tc=round(ns, 3),
        SA=r.sa, QED=r.qed, MW=r.MW, LogP=r.LogP, TPSA=r.TPSA, HBD=r.HBD, HBA=r.HBA, RotB=r.RotB,
        fCsp3=round(rdMolDescriptors.CalcFractionCSP3(mol), 3),
        aromatic_rings=rdMolDescriptors.CalcNumAromaticRings(mol),
        lipinski_violations=lip,
        veber_pass=int(rdMolDescriptors.CalcNumRotatableBonds(mol) <= 10 and Descriptors.TPSA(mol) <= 140),
        n_extra_alerts=len(alerts),
        extra_alerts="; ".join(sorted(set(alerts))[:3]),
        scaffold=r.scaffold, chemotype=r.chemotype
    ))

F = pd.DataFrame(rows)
F['evidence'] = np.where(F.selectivity <= -0.5, 'QSAR + docking + selectivity margin', 'QSAR + docking')

def tier(r):
    if r.n_extra_alerts > 0:
        return 'Tier 3 - structural alert flagged'
    if r.track.startswith('Azole-free') and r.selectivity <= -0.5:
        return 'Tier 1 - lead'
    return 'Tier 2 - supported'

F['tier'] = F.apply(tier, axis=1)
F = F.sort_values(['tier', 'selectivity']).reset_index(drop=True)
F.to_csv("final_candidates.csv", index=False)