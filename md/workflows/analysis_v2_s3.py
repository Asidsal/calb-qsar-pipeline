"""Corrected 5TZ1 CAL-028/CAL-015 analysis for Google Colab (CPU).
Correctly includes 13 AMBER histidine variants omitted by v1 protein selection.
Labels protein residues with original PDB numbering 45–528.
Read-only RunPod S3 access; no GPU required. Requires pre-existing variables s3, bucket.
"""
"""Stream CAL-028/CAL-015 50 ns trajectories from a RunPod S3 volume to Colab.
Requires existing variables s3 and bucket from prior S3 credential cell.
"""
from pathlib import Path
from collections import Counter
import json
import re
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import mdtraj as md

if 's3' not in globals() or 'bucket' not in globals():
    raise RuntimeError('Run the preceding RunPod S3 credentials cell to define s3 and bucket.')

OUTPUT = Path('/content/CAL028_CAL015_results_v2')
OUTPUT.mkdir(parents=True, exist_ok=True)
TEMP = OUTPUT / 'temp_1ns.dcd'

# Paginate: the original first-100 list cannot confirm a completed 50 ns trajectory.
keys = []
for page in s3.get_paginator('list_objects_v2').paginate(Bucket=bucket):
    keys.extend(item['Key'] for item in page.get('Contents', []))

patterns = {
    'CAL028': re.compile(r'^CAL028_5090_results/ns_(\d{3})/trajectory\.dcd$'),
    'CAL015': re.compile(r'^CAL-015/production_50ns_1fs/ns_(\d{3})/trajectory\.dcd$'),
}
selections = {}
for label, rx in patterns.items():
    selections[label] = sorted((int(m.group(1)), key) for key in keys if (m := rx.fullmatch(key)))
    indices = [j for j, _ in selections[label]]
    missing = sorted(set(range(1, 51)) - set(indices))
    print(label, 'DCD segments:', len(indices), 'missing 1–50:', missing)
    if missing:
        raise RuntimeError(f'{label}: missing trajectory segments {missing}; do not call this a complete 50 ns analysis')

choices = {
    'CAL028': ['cal028.prmtop', 'CAL028_5090_results/cal028.prmtop'],
    'CAL015': ['CAL-015/N14_frame_corrected/system.prmtop'],
}
for label in choices:
    candidates = [k for k in choices[label] if k in keys]
    if not candidates:
        print(label, 'available topology candidates:',
              [k for k in keys if k.lower().endswith('.prmtop')][:30])
        raise FileNotFoundError(f'{label}: expected topology not found. Identify correct .prmtop from listing.')
    choices[label] = candidates[0]

summary = {}
for label in ['CAL028', 'CAL015']:
    print('\nAnalyzing', label, 'using', choices[label], flush=True)
    topfile = OUTPUT / f'{label}.prmtop'
    s3.download_file(bucket, choices[label], str(topfile))
    topology = md.load_prmtop(str(topfile), standard_names=False)

    ligand = np.array([a.index for a in topology.atoms if a.residue.name == 'LIG'
                      and a.element and a.element.symbol != 'H'], dtype=int)
    # Identify all polypeptide residues by backbone chemistry, not MDTraj's
    # `residue.is_protein` flag: AMBER's HID/HIE protonation variants may be
    # omitted by that property, losing 13 histidines in these 5TZ1 systems.
    protein_residues = {
        r.index for r in topology.residues
        if r.name not in {'LIG', 'HEM', 'WAT', 'HOH', 'NA', 'CL'}
        and {'N', 'CA', 'C'}.issubset({a.name for a in r.atoms})
    }
    assert len(protein_residues) == 484, (label, len(protein_residues),
        'Expected residues 45–528; check the prepared system.')
    def resid_label(residue):
        if residue.index in protein_residues:
            return f'{residue.name}{residue.index + 45}'  # original PDB residue ID
        return str(residue)
    ca = np.array([a.index for a in topology.atoms
                   if a.name == 'CA' and a.residue.index in protein_residues], dtype=int)
    prot_heavy = np.array([a.index for a in topology.atoms
                          if a.residue.index in protein_residues
                          and a.element and a.element.symbol != 'H'], dtype=int)
    print(f'{label}: recognized {len(protein_residues)} protein residues, '
          f'{len(ca)} C-alpha atoms; mapping protein residue indices +45 to PDB IDs.')
    assert len(ligand) > 5, f'{label} ligand heavy atoms = {len(ligand)}; check residue name'
    assert len(ca) > 300, f'{label} protein C-alpha count unexpectedly small: {len(ca)}'
    fe = next(a.index for a in topology.atoms if a.residue.name == 'HEM' and a.name.upper() == 'FE')
    sg = next(a.index for a in topology.atoms if a.residue.name in {'CYM', 'CYP'} and a.name == 'SG')
    heme_n = {name: next(a.index for a in topology.atoms if a.residue.name == 'HEM' and a.name == name)
              for name in ('NA', 'NB', 'NC', 'ND')}

    ca_sum = np.zeros((len(ca), 3), dtype=np.float64)
    ca_sq = np.zeros((len(ca), 3), dtype=np.float64)
    contacts = Counter()
    hbond_hits = Counter()
    n_hbond_frames = 0
    total_frames = 0
    reference = None
    all_rows = []
    imaging_warnings = set()
    for ns, key in selections[label]:
        s3.download_file(bucket, key, str(TEMP))
        traj = md.load_dcd(str(TEMP), top=topology)
        TEMP.unlink(missing_ok=True)
        assert traj.n_atoms == topology.n_atoms
        try:
            if traj.unitcell_lengths is not None:
                traj.image_molecules(inplace=True)
            else:
                imaging_warnings.add('No unit cell in DCD: PBC not applied')
        except Exception as exc:
            imaging_warnings.add('Imaging failed: ' + str(exc)[:120])

        if reference is None:
            reference = traj[0]
        traj.superpose(reference, frame=0, atom_indices=ca)

        # Critical: ligand RMSD uses protein-frame fit, NOT a ligand-only fit.
        ligand_rmsd = 10 * np.sqrt(np.mean(np.sum((traj.xyz[:, ligand] -
                                               reference.xyz[0, ligand])**2, axis=2), axis=1))
        protein_rmsd = 10 * np.sqrt(np.mean(np.sum((traj.xyz[:, ca] -
                                                reference.xyz[0, ca])**2, axis=2), axis=1))
        ca_xyz = traj.xyz[:, ca, :].astype(np.float64)
        ca_sum += ca_xyz.sum(axis=0)
        ca_sq += np.square(ca_xyz).sum(axis=0)
        center = ca_xyz.mean(axis=1, keepdims=True)
        ca_rg = 10 * np.sqrt(np.mean(np.sum((ca_xyz-center)**2, axis=2), axis=1))
        fe_s = 10 * np.linalg.norm(traj.xyz[:, fe]-traj.xyz[:, sg], axis=1)
        fe_n = {nm: 10*np.linalg.norm(traj.xyz[:, fe]-traj.xyz[:, j], axis=1)
                for nm, j in heme_n.items()}

        neighbors = md.compute_neighbors(traj, 0.45, query_indices=ligand,
                                         haystack_indices=prot_heavy, periodic=False)
        for frame_neigh in neighbors:
            res = {resid_label(topology.atom(int(idx)).residue) for idx in frame_neigh}
            contacts.update(res)
        # Protein-ligand hydrogen bonds, excluding water; sample 1 frame in 5.
        try:
            sampled = traj[::5]
            detected = md.wernet_nilsson(sampled, exclude_water=True, periodic=False)
            for frame_triplets in detected:
                hits = set()
                for d, h, a in frame_triplets:
                    donor, acceptor = topology.atom(int(d)), topology.atom(int(a))
                    if ((donor.residue.name == 'LIG' and acceptor.residue.index in protein_residues) or
                        (acceptor.residue.name == 'LIG' and donor.residue.index in protein_residues)):
                        hits.add(f'{resid_label(donor.residue)}:{donor.name} -> '
                                 f'{resid_label(acceptor.residue)}:{acceptor.name}')
                hbond_hits.update(hits)
            n_hbond_frames += len(detected)
        except Exception as exc:
            imaging_warnings.add('Hydrogen-bond analysis skipped: ' + str(exc)[:120])

        n = traj.n_frames
        # DCD files are 1 ns segments; assign uniformly spaced times within each.
        times_ns = (ns-1) + (np.arange(1, n+1) / n)
        for j in range(n):
            all_rows.append({'time_ns': float(times_ns[j]),
                             'protein_CA_RMSD_A': float(protein_rmsd[j]),
                             'ligand_protein_aligned_RMSD_A': float(ligand_rmsd[j]),
                             'protein_CA_Rg_A': float(ca_rg[j]),
                             'Fe_S_A': float(fe_s[j]),
                             **{f'Fe_{nm}_A': float(v[j]) for nm, v in fe_n.items()}})
        total_frames += n
        if ns % 5 == 0:
            print(f'{label}: {ns}/50 ns, frames={total_frames}', flush=True)

    df = pd.DataFrame(all_rows)
    df.to_csv(OUTPUT / f'{label}_time_series.csv', index=False)
    rmsf_A = 10 * np.sqrt(np.maximum(0, np.sum(ca_sq/total_frames -
                                              np.square(ca_sum/total_frames), axis=1)))
    pd.DataFrame([{'residue': resid_label(topology.atom(int(j)).residue),
                   'CA_RMSF_A': float(val)} for j, val in zip(ca, rmsf_A)]).to_csv(
                       OUTPUT / f'{label}_CA_RMSF.csv', index=False)
    pd.DataFrame([{'residue': key, 'occupancy_percent': 100*count/total_frames}
                  for key, count in contacts.most_common()]).to_csv(
                      OUTPUT / f'{label}_ligand_contacts_4p5A.csv', index=False)
    pd.DataFrame([{'hydrogen_bond': key, 'sampled_occupancy_percent': 100*count/max(1,n_hbond_frames)}
                  for key, count in hbond_hits.most_common()]).to_csv(
                      OUTPUT / f'{label}_ligand_hbonds_sampled.csv', index=False)

    summary[label] = {'frames': total_frames, 'protein_CA_atoms': len(ca),
                      'residue_labeling': 'PDB numbering; prepared protein residue index + 45',
                      'mean_ligand_RMSD_A': float(df['ligand_protein_aligned_RMSD_A'].mean()),
                      'mean_protein_CA_RMSD_A': float(df['protein_CA_RMSD_A'].mean()),
                      'mean_Fe_S_A': float(df['Fe_S_A'].mean()),
                      'imaging_warnings': sorted(imaging_warnings),
                      'source_topology': choices[label],
                      'top_contacts': [{'residue': k, 'occupancy_percent': 100*v/total_frames}
                                       for k,v in contacts.most_common(10)]}

(OUTPUT / 'summary.json').write_text(json.dumps(summary, indent=2))
fig, axes = plt.subplots(3, 1, figsize=(10, 11), sharex=True)
for label in ('CAL028', 'CAL015'):
    df = pd.read_csv(OUTPUT / f'{label}_time_series.csv')
    for ax, col in zip(axes,
        ['ligand_protein_aligned_RMSD_A','protein_CA_RMSD_A','Fe_S_A']):
        ax.plot(df['time_ns'], df[col], label=label, alpha=0.75, linewidth=0.8)
for ax, ylabel in zip(axes, ['Ligand RMSD (Å), protein-aligned',
                             'Protein CA RMSD (Å)', 'Fe–S distance (Å)']):
    ax.set_ylabel(ylabel)
    ax.legend()
axes[-1].set_xlabel('Production time (ns)')
fig.tight_layout()
fig.savefig(OUTPUT / 'CAL028_CAL015_trajectory_comparison.png', dpi=220)
plt.show()

fig2, ax2 = plt.subplots(figsize=(10,4))
for label in ('CAL028','CAL015'):
    df2 = pd.read_csv(OUTPUT / f'{label}_CA_RMSF.csv')
    ax2.plot(range(1,len(df2)+1), df2['CA_RMSF_A'], label=label, linewidth=0.9)
ax2.set(xlabel='Protein C-alpha ordinal position (residues 45-528)', ylabel='CA RMSF (Å)')
ax2.legend()
fig2.tight_layout()
fig2.savefig(OUTPUT / 'CAL028_CAL015_protein_CA_RMSF.png', dpi=220)
plt.show()

fig3, axs3 = plt.subplots(1, 2, figsize=(12,5))
for ax3, label in zip(axs3, ('CAL028','CAL015')):
    d = pd.read_csv(OUTPUT / f'{label}_ligand_contacts_4p5A.csv').head(12)
    ax3.barh(d['residue'][::-1], d['occupancy_percent'][::-1])
    ax3.set(title=label, xlabel='Contact occupancy (%)', xlim=(0,100))
fig3.tight_layout()
fig3.savefig(OUTPUT / 'CAL028_CAL015_contact_occupancy.png', dpi=220)
plt.show()
print('RESULTS:', OUTPUT)
print(json.dumps(summary, indent=2)[:5000])

# Small analysis outputs only: do not package full trajectories or credentials.
import zipfile
zip_path = OUTPUT.parent / 'CAL028_CAL015_analysis_results_v2.zip'
with zipfile.ZipFile(zip_path, 'w', compression=zipfile.ZIP_DEFLATED) as zf:
    for f in sorted(OUTPUT.iterdir()):
        if f.is_file() and f.suffix in {'.csv','.png','.json'}:
            zf.write(f, arcname=f.name)
print('RESULTS ZIP:', zip_path)
