load CAL015_median_pose_protein_heme_ligand.pdb, complex
hide everything, all
show cartoon, polymer.protein
color gray70, polymer.protein
show sticks, resn LIG or resn HEM
color yellow, resn LIG
color salmon, resn HEM
select nearby, byres (polymer.protein within 4.5 of resn LIG)
show sticks, nearby
color cyan, nearby
bg_color white
orient resn LIG
set ray_opaque_background, 0
png CAL015_binding_site_3D.png, width=1800, height=1400, dpi=300, ray=1
