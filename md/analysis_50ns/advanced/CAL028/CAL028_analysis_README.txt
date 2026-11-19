CAL028: 50 ns standalone trajectory analysis
The analysis uses 5000 frames across 50 1-ns DCD segments; ligand pose RMSD is protein-CA aligned against (i) first production frame and (ii) original SDF/PDB docking pose.
SASA uses heavy-atom Shrake-Rupley radii with 80 sphere points, 4 snapshots/ns. Ligand isolated SASA uses the SAME snapshot conformations, removing the receptor; ligand buried SASA = isolated minus bound. Explicit water is omitted in SASA surfaces.
Contacts use 4.5 A ligand-protein heavy atom proximity over all 5000 frames. C/S packing proxy is a geometric heuristic, not an energy estimate.
Wernet-Nilsson direct ligand-protein H-bonds use every fifth frame, excluding water.
SDF stereo assignments are inferred from original three-dimensional reference. Alkene dihedral traces are from the corresponding topological atoms in the MD coordinates.
Representative PDB/SDF files are protein-aligned simulation snapshots. Use a molecular viewer (PyMOL/ChimeraX) for publication 3D poses.
All figures describe one trajectory per complex, not replicate variability.
