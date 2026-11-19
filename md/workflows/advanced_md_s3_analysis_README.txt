CAL-028 / CAL-015 Advanced Molecular Dynamics Analysis
===================================================

Open CAL028_CAL015_Advanced_MD_Analysis.ipynb in Google Colab.
Select CPU runtime and Run all. The notebook has the original docking references
embedded, connects securely to the RunPod S3 volume (asks you for S3 credentials
using hidden input), downloads 50 DCD segments each compound one at a time, and
saves publication figures and CSVs.

INPUT
-----
RunPod bucket: lf87cdnmsz (volume cal028-md, EU-RO-1)
CAL028 DCD: CAL028_5090_results/ns_001...ns_050/trajectory.dcd
CAL028 prmtop: cal028.prmtop
CAL015 DCD: CAL-015/production_50ns_1fs/ns_001...ns_050/trajectory.dcd
CAL015 prmtop: CAL-015/N14_frame_corrected/system.prmtop

ANALYSES
--------
1. SASA (bound and free ligand, protein, polar heavy atom contributions, burial)
2. Direct H-bond count and timelines; 4.5-A contact and packing proxy heatmaps
3. Protein-aligned RMSD to first production frame AND original docked pose
4. Ligand-heme Fe nearest atom/N distance, Fe-S bond, CA RMSF
5. CAL-015 substituted alkene E/Z geometric dihedral monitoring
6. First/median/final protein-heme-ligand PDB and mapped ligand SDF

Standalone figures already calculated from the previous v2 numerical output:
 individual_v2_figures/CAL028_50ns_Individual_MD_Summary.png
 individual_v2_figures/CAL015_50ns_Individual_MD_Summary.png

New SASA/docking-reference/time-resolved figures require Colab to read the
actual DCD trajectories. They are not in this ZIP until the notebook runs.

This notebook uses a CPU only; it does NOT launch MD, alter RunPod data or
restart a paid GPU. Do not upload your S3 credentials here.
