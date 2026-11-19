"""Individual advanced MD analyses for CAL-028 and CAL-015, 50 ns each.

Run in Colab CPU. Input: read-only S3 objects on private RunPod network volume.
Contains no RunPod secrets. References are original docked PDB/SDF files.
Scientific modules: SASA, time-resolved polar/nonpolar interactions,
original docking pose RMSD, E/Z monitoring, representative PDB/SDF exports.

The code computes NEW results only when executed with actual DCDs.
"""
import base64
import csv
import gzip
import json
import math
import os
import re
import sys
import zipfile
from pathlib import Path
from collections import Counter, defaultdict

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import mdtraj as md
import networkx as nx
from rdkit import Chem
from rdkit.Chem import rdMolTransforms

S3_BUCKET = 'lf87cdnmsz'
S3_ENDPOINT = 'https://s3api-eu-ro-1.runpod.io'
REGION = 'eu-ro-1'
OUTPUT = Path('/content/CAL028_CAL015_advanced_results') if Path('/content').exists() else Path.cwd() / 'CAL028_CAL015_advanced_results'
REF_ROOT = Path('/content/CAL_MD_reference_poses') if Path('/content').exists() else Path.cwd() / 'CAL_MD_reference_poses'
OUTPUT.mkdir(parents=True, exist_ok=True)
TEMP = OUTPUT / '_temp_segment.dcd'
# 4 SASA frames per ns; Wernet–Nilsson hydrogen bonds on 20 frames/ns.
SASA_FRAME_STRIDE = 25
HBOND_FRAME_STRIDE = 5
SASA_SPHERE_POINTS = 80
CONTACT_CUTOFF_NM = 0.45
EXPECTED_SEGMENTS = 50
LABELS = ['CAL028', 'CAL015']
S3_TOPOLOGY = {
    'CAL028': ['cal028.prmtop', 'CAL028_5090_results/cal028.prmtop'],
    'CAL015': ['CAL-015/N14_frame_corrected/system.prmtop'],
}
S3_SEGMENT_REGEX = {
    'CAL028': re.compile(r'^CAL028_5090_results/ns_(\d{3})/trajectory\.dcd$'),
    'CAL015': re.compile(r'^CAL-015/production_50ns_1fs/ns_(\d{3})/trajectory\.dcd$'),
}
SOURCE_LABEL = {'CAL028': 'CAL-028', 'CAL015': 'CAL-015'}


def get_s3_client():
    # Reuse the secure client from a preceding Colab cell when available.
    existing = globals().get('s3')
    if existing is not None:
        return existing, globals().get('bucket', S3_BUCKET)
    import boto3
    from getpass import getpass
    from botocore.config import Config
    uid = getpass('RunPod S3 User ID (hidden): ')
    secret = getpass('RunPod S3 secret (hidden): ')
    # Do not persist or print credentials anywhere.
    client = boto3.client('s3', endpoint_url=S3_ENDPOINT, region_name=REGION,
                          aws_access_key_id=uid, aws_secret_access_key=secret,
                          config=Config(signature_version='s3v4',
                                        s3={'addressing_style':'path'}))
    del uid, secret
    return client, S3_BUCKET


def list_s3_objects(client, bucket_name):
    key_list=[]
    for page in client.get_paginator('list_objects_v2').paginate(Bucket=bucket_name):
        key_list.extend(item['Key'] for item in page.get('Contents', []))
    key_set=set(key_list)
    choices={}
    for label in LABELS:
        found=sorted((int(m.group(1)), key)
                     for key in key_list
                     if (m:=S3_SEGMENT_REGEX[label].fullmatch(key)))
        ids=[n for n,_ in found]
        missing=sorted(set(range(1,EXPECTED_SEGMENTS+1))-set(ids))
        if missing:
            raise RuntimeError(f'{label}: missing 1ns trajectory segments {missing}; no partial 50ns result will be reported')
        extra=sorted(set(ids)-set(range(1,EXPECTED_SEGMENTS+1)))
        if extra: print(f'Note: {label} also has additional trajectory segments {extra}; analyzing 1–50 ns only')
        choices[label]={'segments':[(n,k) for n,k in found if 1<=n<=EXPECTED_SEGMENTS]}
        top_candidates=[t for t in S3_TOPOLOGY[label] if t in key_set]
        if not top_candidates:
            raise FileNotFoundError(f'{label}: missing expected AMBER prmtop; known candidates {S3_TOPOLOGY[label]}')
        choices[label]['topology']=top_candidates[0]
        print(f'{label}: all 50 DCD segments found; topology {top_candidates[0]}')
    return choices


def original_ca_and_ligand(reference_pdb):
    cas={}
    for line in Path(reference_pdb).read_text(errors='replace').splitlines():
        if line.startswith('ATOM  ') and line[12:16].strip()=='CA':
            resid=int(line[22:26]); chain=line[21].strip()
            if chain=='A' and resid not in cas:
                cas[resid]=np.array([float(line[i:i+8]) for i in (30,38,46)],dtype=float)/10
    return cas


def kabsch_to_frame(P,Q):
    """Return rigid transform aligning rows P to rows Q, both in nm."""
    cp=P.mean(axis=0);cq=Q.mean(axis=0)
    U,s,Vt=np.linalg.svd((P-cp).T@(Q-cq))
    d=np.sign(np.linalg.det(U@Vt))
    R=U@np.diag([1.0,1.0,d])@Vt
    return R, cq-cp@R


def identify_system(top):
    protein={r.index for r in top.residues
             if r.name not in {'HEM','LIG','WAT','HOH','NA','CL','Na+','Cl-'}
             and {'N','CA','C'}.issubset({a.name for a in r.atoms})}
    ca=np.array([a.index for a in top.atoms if a.name=='CA' and a.residue.index in protein],dtype=int)
    if len(protein)!=484 or len(ca)!=484:
        raise AssertionError(f'Expected 484 protein residues / CA atoms, got {len(protein)}/{len(ca)}; inspect prepared topology')
    protein_heavy=np.array([a.index for a in top.atoms
                            if a.residue.index in protein and a.element and a.element.symbol!='H'],dtype=int)
    heme_heavy=np.array([a.index for a in top.atoms
                          if a.residue.name=='HEM' and a.element and a.element.symbol!='H'], dtype=int)
    lig_heavy=np.array([a.index for a in top.atoms
                         if a.residue.name=='LIG' and a.element and a.element.symbol!='H'],dtype=int)
    if len(lig_heavy)!=32:
        raise AssertionError(f'Expected 32 ligand heavy atoms, got {len(lig_heavy)}')
    lig_nonpolar=np.array([i for i in lig_heavy if top.atom(int(i)).element.symbol in {'C','S','Cl','Br','I','F'}],dtype=int)
    prot_nonpolar=np.array([i for i in protein_heavy if top.atom(int(i)).element.symbol in {'C','S'}],dtype=int)
    fe=next(a.index for a in top.atoms if a.residue.name=='HEM' and a.name.upper()=='FE')
    sg=next(a.index for a in top.atoms if a.name=='SG' and a.residue.name in {'CYM','CYP'})
    heme_n={nm: next(a.index for a in top.atoms if a.residue.name=='HEM' and a.name==nm)
            for nm in ['NA','NB','NC','ND']}
    return locals()


def pdb_res_label(residue, protein_index_set):
    if residue.index in protein_index_set:
        return f'{residue.name}{residue.index+45}'
    return f'{residue.name}{residue.resSeq}'


def ligand_graph_mapping(topology, original_mol, topology_lig_indices, original_xyz_aligned, first_xyz):
    """Match ligand covalent heavy-atom graphs; choose symmetry mapping closest to original pose."""
    gref=nx.Graph()
    for atom in original_mol.GetAtoms():
        if atom.GetSymbol()!='H':
            gref.add_node(atom.GetIdx(), element=atom.GetSymbol().upper())
    for bond in original_mol.GetBonds():
        i,j=bond.GetBeginAtomIdx(),bond.GetEndAtomIdx()
        if i in gref and j in gref: gref.add_edge(i,j)
    gt=nx.Graph()
    ligset=set(map(int,topology_lig_indices))
    for i in topology_lig_indices:
        a=topology.atom(int(i))
        gt.add_node(int(i), element=a.element.symbol.upper())
    for a,b in topology.bonds:
        if a.index in ligset and b.index in ligset:gt.add_edge(a.index,b.index)
    matcher=nx.algorithms.isomorphism.GraphMatcher(
        gref,gt,node_match=lambda a,b:a['element']==b['element'])
    best=None; bestscore=float('inf'); tested=0
    refids=sorted(gref.nodes)
    for mapping in matcher.isomorphisms_iter():
        dst=np.array([mapping[i] for i in refids],dtype=int)
        score=float(np.sqrt(np.mean(np.sum((first_xyz[dst]-original_xyz_aligned[refids])**2,axis=1))))
        if score<bestscore: best=dict(mapping);bestscore=score
        tested+=1
        if tested>=512:break
    if best is None:
        raise ValueError('Original SDF and prepared topology have different ligand heavy-atom connectivity; docking pose RMSD cannot be calculated safely.')
    print(f'    Ligand covalent graph matches: {tested}, best first-frame distance RMSD {bestscore*10:.3f} Å')
    return best


def load_reference(label, top, info, frame):
    stem=SOURCE_LABEL[label]
    ref_pdb=REF_ROOT/f'{stem}_complex.pdb'
    ref_sdf=REF_ROOT/f'{stem}_ligand_pose.sdf'
    if not ref_pdb.exists() or not ref_sdf.exists():
        raise FileNotFoundError(f'Missing original pose references: {ref_pdb}, {ref_sdf}')
    ca_ref=original_ca_and_ligand(ref_pdb)
    expected_ids=[int(top.atom(int(i)).residue.index+45) for i in info['ca']]
    assert all(r in ca_ref for r in expected_ids), 'Original PDB does not contain all 484 matching CA residue IDs'
    P=np.array([ca_ref[i] for i in expected_ids])
    Q=frame.xyz[0,info['ca'],:]
    R,t=kabsch_to_frame(P,Q)
    protein_ca_fit_rms=float(np.sqrt(np.mean(np.sum(((P@R+t)-Q)**2,axis=1)))*10)
    original_mol=Chem.MolFromMolFile(str(ref_sdf),removeHs=False)
    if original_mol is None: raise ValueError(f'Invalid original SDF: {ref_sdf}')
    conformer=original_mol.GetConformer()
    xyz=np.array(conformer.GetPositions(),dtype=float)/10
    xyzaligned=xyz@R+t
    mapping=ligand_graph_mapping(top, original_mol, info['lig_heavy'], xyzaligned,frame.xyz[0])
    rdids=np.array(sorted(mapping),dtype=int)
    mapped_ids=np.array([mapping[i] for i in rdids],dtype=int)
    Chem.AssignStereochemistryFrom3D(original_mol)
    stereo_bonds=[]
    for b in original_mol.GetBonds():
        if b.GetBondType()!=Chem.BondType.DOUBLE or b.GetStereo() not in (Chem.BondStereo.STEREOE,Chem.BondStereo.STEREOZ):continue
        i,j=b.GetBeginAtomIdx(),b.GetEndAtomIdx()
        st=list(b.GetStereoAtoms())
        if len(st)==2 and all(q in mapping for q in [st[0],i,j,st[1]]):
            n_i={a.GetIdx() for a in original_mol.GetAtomWithIdx(i).GetNeighbors()}
            n_j={a.GetIdx() for a in original_mol.GetAtomWithIdx(j).GetNeighbors()}
            # RDKit's stereo-atom list may be in reverse bond order after
            # AssignStereochemistryFrom3D. Establish a connected torsion.
            if st[0] in n_i and st[1] in n_j:
                quadruplet=[st[0],i,j,st[1]]
            elif st[0] in n_j and st[1] in n_i:
                quadruplet=[st[0],j,i,st[1]]
            else:
                raise ValueError(f'RDKit stereo atoms {st} do not neighbor C=C {i}-{j}')
            orig_dihedral=rdMolTransforms.GetDihedralDeg(conformer,*quadruplet)
            stereo_bonds.append({'bond_atoms_SDF':[i,j],
                                 'reference_geometry':str(b.GetStereo()),
                                 'quadruplet_topology':[mapping[k] for k in quadruplet],
                                 'original_dihedral_deg':float(orig_dihedral)})
    return {'mol':original_mol,'docked_xyz':xyzaligned,'reference_ids':rdids,
            'mapped_ids':mapped_ids,'atom_map':mapping,'stereo_bonds':stereo_bonds,
            'protein_ca_fit_A':protein_ca_fit_rms,'docking_frame_R':R.tolist(),'docking_frame_t_nm':t.tolist()}


def time_by_segment(ns, n_frames):
    return (ns-1)+np.arange(1,n_frames+1,dtype=float)/n_frames


def _plot_curve(ax,x,y,title,ylabel):
    ax.plot(x,y,linewidth=0.75)
    ax.set_title(title)
    ax.set_xlabel('Time (ns)')
    ax.set_ylabel(ylabel)
    ax.grid(alpha=.2)


def _save_png(fig,filename):
    fig.tight_layout()
    fig.savefig(filename,dpi=280,bbox_inches='tight')
    plt.close(fig)


def save_snapshot(label,tag,segment_key,frame_index,client,bucket_name,top,info,original_mol,mapping,reference,odir):
    print(f'    Exporting {tag} representative: segment {segment_key}, frame {frame_index}',flush=True)
    client.download_file(bucket_name,segment_key,str(TEMP))
    segment=md.load_dcd(str(TEMP),top=top)
    TEMP.unlink(missing_ok=True)
    segment.superpose(reference,atom_indices=info['ca'])
    snap=segment[frame_index]
    subset_idx=np.unique(np.concatenate([info['protein_heavy'], info['heme_heavy'],
                                         info['lig_heavy']]))
    outp=odir/f'{label}_{tag}_protein_heme_ligand.pdb'
    snap.atom_slice(subset_idx).save_pdb(str(outp))
    # AMBER topologies may store sequential protein residue IDs (1–484).
    # Remap those to original 5TZ1 IDs 45–528 in exported PDBs.
    lines=outp.read_text().splitlines()
    ca_lines=[x for x in lines if x.startswith(('ATOM  ','HETATM'))
              and x[12:16].strip()=='CA' and x[17:20].strip() not in ('LIG','HEM')]
    if ca_lines and int(ca_lines[0][22:26])==1:
        corrected=[]
        for line in lines:
            if line.startswith(('ATOM  ','HETATM')) and line[17:20].strip() not in ('LIG','HEM'):
                old=int(line[22:26]);line=line[:22]+f'{old+44:4d}'+line[26:]
            corrected.append(line)
        outp.write_text('\n'.join(corrected)+'\n')
    # Ligand SDF retains the original bond order and atom labels, but uses the observed MD heavy-atom geometry.
    mol=Chem.RemoveHs(original_mol)
    conf=mol.GetConformer()
    for origidx,topidx in mapping.items():
        if origidx>=mol.GetNumAtoms():
            # RDKit usually stores all heavy atoms before hydrogens in the original SDF.
            continue
        p=snap.xyz[0,topidx]*10
        conf.SetAtomPosition(origidx,(float(p[0]),float(p[1]),float(p[2])))
    sdf=odir/f'{label}_{tag}_ligand_heavy.sdf'
    wr=Chem.SDWriter(str(sdf));wr.write(mol);wr.close()
    return {'pdb':str(outp.name),'sdf':str(sdf.name)}


def run_one(label,selection,client,bucket_name):
    print(f'\n========== {label}: independent 50ns analysis ==========',flush=True)
    odir=OUTPUT/label
    odir.mkdir(exist_ok=True)
    top_path=odir/f'{label}_system.prmtop'
    client.download_file(bucket_name,selection['topology'],str(top_path))
    top=md.load_prmtop(str(top_path),standard_names=False)
    info=identify_system(top)
    # Avoid storing the full AMBER topology (10+ MB) in the output ZIP; original remains on RunPod S3.
    top_path.unlink(missing_ok=True)
    meta={'system':label,'topology_source':selection['topology'],
          'trajectory_files':[key for _,key in selection['segments']],
          'analysis_settings':{'sasa_stride_frames':SASA_FRAME_STRIDE,
             'sasa_sphere_points':SASA_SPHERE_POINTS,'hbond_stride_frames':HBOND_FRAME_STRIDE,
             'contact_cutoff_A':CONTACT_CUTOFF_NM*10,
             'sasa_waters_excluded':True},'protein_CA_count':len(info['ca'])}
    # Precompute selection mapping for mdtraj SASA on heavy atom subsets.
    system_sasa=np.sort(np.concatenate([info['protein_heavy'],info['heme_heavy'],info['lig_heavy']]))
    sys_pos={int(x):p for p,x in enumerate(system_sasa)}
    lig_pos=[sys_pos[int(x)] for x in info['lig_heavy']]
    prot_pos=[sys_pos[int(x)] for x in info['protein_heavy']]
    polar_pos=[p for p,i in zip(lig_pos,info['lig_heavy']) if top.atom(int(i)).element.symbol in {'N','O','S'}]
    other_pos=[p for p in lig_pos if p not in set(polar_pos)]
    frame_records=[];sasa_records=[];contact_events=[];hydrogen_events=[];hbond_count_records=[]
    per_frame_contacts=[]; ca_sum=np.zeros((len(info['ca']),3),dtype=np.float64)
    ca_sq=np.zeros_like(ca_sum)
    reference=None; docking=None; N=0; seg_lookup={}
    total_hb_sampled=0
    for ns,key in selection['segments']:
        client.download_file(bucket_name,key,str(TEMP))
        traj=md.load_dcd(str(TEMP),top=top)
        TEMP.unlink(missing_ok=True)
        assert traj.n_atoms==top.n_atoms, f'{label} segment {ns}: topology count mismatch'
        assert traj.n_frames==100, f'{label} segment {ns}: expected 100 DCD frames, got {traj.n_frames}'
        if reference is None:
            reference=traj[0]
            docking=load_reference(label,top,info,reference)
            meta['original_pose_CA_alignment_RMSD_A']=docking['protein_ca_fit_A']
            meta['original_pose_mapping']={str(k):int(v) for k,v in docking['atom_map'].items()}
            meta['stereochemistry_from_original_3D_SDF']=[{k:v for k,v in d.items() if k!='quadruplet_topology'} for d in docking['stereo_bonds']]
            print(f'    Original protein CA fit RMSD = {docking["protein_ca_fit_A"]:.3f} Å')
            if docking['stereo_bonds']:print('    Observed SDF double bond stereochemistry:',meta['stereochemistry_from_original_3D_SDF'])
        try:
            if traj.unitcell_lengths is not None:traj.image_molecules(inplace=True)
            else:print('WARN: DCD unit cell absent: using non-PBC measurements')
        except Exception as exc:
            print('WARN: image_molecules failed:',str(exc)[:110])
        traj.superpose(reference,atom_indices=info['ca'])
        times=time_by_segment(ns,traj.n_frames)
        seg_lookup[ns]=(key,traj.n_frames)
        coord=traj.xyz
        ligand_rms0=10*np.sqrt(np.mean(np.sum(
            (coord[:,info['lig_heavy']]-reference.xyz[0,info['lig_heavy']])**2,axis=2),axis=1))
        ligand_rms_dock=10*np.sqrt(np.mean(np.sum(
            (coord[:,docking['mapped_ids']]-docking['docked_xyz'][docking['reference_ids']])**2,axis=2),axis=1))
        ca_xyz=coord[:,info['ca']].astype(np.float64)
        protein_rms=10*np.sqrt(np.mean(np.sum((ca_xyz-reference.xyz[0,info['ca']])**2,axis=2),axis=1))
        ca_sum+=ca_xyz.sum(axis=0)
        ca_sq+=np.square(ca_xyz).sum(axis=0)
        Fe_S_A=10*np.linalg.norm(coord[:,info['fe']]-coord[:,info['sg']],axis=1)
        Fe_min_lig_A=10*np.linalg.norm(coord[:,info['lig_heavy']]-coord[:,info['fe'],None],axis=2).min(axis=1)
        ligand_N=[a for a in info['lig_heavy'] if top.atom(int(a)).element.symbol=='N']
        Fe_min_N_A=(10*np.linalg.norm(coord[:,ligand_N]-coord[:,info['fe'],None],axis=2).min(axis=1)
                    if ligand_N else np.full(traj.n_frames,np.nan))
        stereo_series=[]
        for b in docking['stereo_bonds']:
            vals=np.degrees(md.compute_dihedrals(traj,np.array([b['quadruplet_topology']]),periodic=False).flatten())
            stereo_series.append(vals)
        # Radius-of-gyration of C-alpha atoms (same convention as the previous analysis).
        center=ca_xyz.mean(axis=1,keepdims=True)
        CaRgA=10*np.sqrt(np.mean(np.sum((ca_xyz-center)**2,axis=2),axis=1))
        for fi in range(traj.n_frames):
            row={'time_ns':float(times[fi]),'ns_segment':ns,'frame_in_segment':fi,
                 'ligand_RMSD_vs_production_A':float(ligand_rms0[fi]),
                 'ligand_RMSD_vs_original_docking_A':float(ligand_rms_dock[fi]),
                 'protein_CA_RMSD_A':float(protein_rms[fi]),
                 'protein_CA_Rg_A':float(CaRgA[fi]),
                 'Fe_S_A':float(Fe_S_A[fi]),
                 'Fe_ligand_min_heavy_A':float(Fe_min_lig_A[fi]),
                 'Fe_ligand_min_N_A':float(Fe_min_N_A[fi])}
            for bi,v in enumerate(stereo_series):row[f'doublebond_{bi+1}_dihedral_deg']=float(v[fi])
            frame_records.append(row)
        # Residue contacts, 100% of trajectory frames, 4.5 Å cutoff.
        contact_per_frame=md.compute_neighbors(traj,CONTACT_CUTOFF_NM,
                  query_indices=info['lig_heavy'],haystack_indices=info['protein_heavy'],periodic=False)
        hydro_contacts_per_frame=md.compute_neighbors(traj,CONTACT_CUTOFF_NM,
                  query_indices=info['lig_nonpolar'],haystack_indices=info['prot_nonpolar'],periodic=False)
        for fi,(cset,hset) in enumerate(zip(contact_per_frame,hydro_contacts_per_frame)):
            allcontact={pdb_res_label(top.atom(int(i)).residue,info['protein']) for i in cset}
            hydrocontact={pdb_res_label(top.atom(int(i)).residue,info['protein']) for i in hset}
            for r in allcontact:contact_events.append((float(times[fi]),ns,r,'all_heavy'))
            for r in hydrocontact:contact_events.append((float(times[fi]),ns,r,'C_S_packing_proxy'))
        # H-bonds: Wernet–Nilsson geometric criteria, exclude waters, every fifth frame.
        hbframes=traj[::HBOND_FRAME_STRIDE]
        hbstep_times=times[::HBOND_FRAME_STRIDE]
        hbset=md.wernet_nilsson(hbframes,exclude_water=True,periodic=False)
        for fi,triplets in enumerate(hbset):
            hits=set()
            for donor,H,acceptor in triplets:
                d=top.atom(int(donor));a=top.atom(int(acceptor))
                if ((d.residue.name=='LIG' and a.residue.index in info['protein'])
                    or (a.residue.name=='LIG' and d.residue.index in info['protein'])):
                    donorname=f'{pdb_res_label(d.residue,info["protein"])}:{d.name}'
                    accname=f'{pdb_res_label(a.residue,info["protein"])}:{a.name}'
                    hits.add(f'{donorname} -> {accname}')
            for hit in hits:hydrogen_events.append((float(hbstep_times[fi]),ns,hit))
            hbond_count_records.append((float(hbstep_times[fi]),ns,len(hits)))
        total_hb_sampled+=len(hbset)
        # Sample SASA on 4 frames/ns; no explicit waters in SASA surface geometry.
        sample_ids=np.arange(0,traj.n_frames,SASA_FRAME_STRIDE)
        subframes=traj[sample_ids]
        complex_sub=subframes.atom_slice(system_sasa)
        free_sub=subframes.atom_slice(info['lig_heavy'])
        # Å²: 1 nm² = 100 Å². Sampled atom radii based on MDTraj's Shrake–Rupley defaults.
        sasa_nm2=md.shrake_rupley(complex_sub,n_sphere_points=SASA_SPHERE_POINTS,mode='atom')
        ligand_free_nm2=md.shrake_rupley(free_sub,n_sphere_points=SASA_SPHERE_POINTS,mode='atom')
        for j,idx in enumerate(sample_ids):
            bound=float(sasa_nm2[j,lig_pos].sum()*100)
            free=float(ligand_free_nm2[j].sum()*100)
            sasa_records.append({'time_ns':float(times[idx]),
                 'protein_SASA_in_complex_A2':float(sasa_nm2[j,prot_pos].sum()*100),
                 'ligand_bound_SASA_A2':bound,'ligand_isolated_SASA_A2':free,
                 'ligand_buried_SASA_A2':max(0.,free-bound),
                 'ligand_bound_NOS_SASA_A2':float(sasa_nm2[j,polar_pos].sum()*100),
                 'ligand_bound_other_heavy_SASA_A2':float(sasa_nm2[j,other_pos].sum()*100)})
        N+=traj.n_frames
        if ns%5==0:print(f'    {label}: processed {ns}/50 ns, {N} frames',flush=True)
        del traj
    assert N==5000, f'{label}: expected 5000 total frames, got {N}'
    df=pd.DataFrame(frame_records);ds=pd.DataFrame(sasa_records)
    df.to_csv(odir/f'{label}_pose_dynamics_timeseries.csv',index=False)
    ds.to_csv(odir/f'{label}_SASA_timeseries.csv',index=False)
    # Recalculate CA RMSF; full 484 residue mapping.
    rmsf_A=10*np.sqrt(np.maximum(0.,np.sum(ca_sq/N-np.square(ca_sum/N),axis=1)))
    rmsf=pd.DataFrame({'pdb_resid':[int(top.atom(int(i)).residue.index)+45 for i in info['ca']],
                      'residue':[pdb_res_label(top.atom(int(i)).residue,info['protein']) for i in info['ca']],
                      'CA_RMSF_A':rmsf_A})
    rmsf.to_csv(odir/f'{label}_protein_CA_RMSF.csv',index=False)
    # Contact and HB occupancy use frame presence, not number of underlying atom pairs.
    ce=pd.DataFrame(contact_events,columns=['time_ns','ns_segment','residue','contact_type'])
    ce.to_csv(odir/f'{label}_contact_frame_events.csv',index=False)
    contact_occup=(ce.groupby(['contact_type','residue']).size().rename('n_frames').reset_index())
    contact_occup['occupancy_percent']=100*contact_occup.n_frames/N
    contact_occup.sort_values(['contact_type','occupancy_percent'],ascending=[True,False],inplace=True)
    contact_occup.to_csv(odir/f'{label}_contacts_occupancy.csv',index=False)
    contact_timeline=(ce.groupby(['contact_type','residue','ns_segment']).size()
                       .rename('n_frames').reset_index())
    contact_timeline['occupancy_percent']=100*contact_timeline.n_frames/100
    contact_timeline.to_csv(odir/f'{label}_contact_1ns_timeline.csv',index=False)
    hb=pd.DataFrame(hydrogen_events,columns=['time_ns','ns_segment','hydrogen_bond'])
    hb.to_csv(odir/f'{label}_Hbond_sampled_events.csv',index=False)
    hbc=pd.DataFrame(hbond_count_records,columns=['time_ns','ns_segment','n_hbonds'])
    hbc.to_csv(odir/f'{label}_Hbond_count_timeseries.csv',index=False)
    if not hb.empty:
        hbocc=(hb.groupby('hydrogen_bond').size().rename('n_sampled_frames').reset_index())
        hbocc['occupancy_percent']=100*hbocc.n_sampled_frames/total_hb_sampled
        hbocc.sort_values('occupancy_percent',ascending=False,inplace=True)
        hbtime=hb.groupby(['hydrogen_bond','ns_segment']).size().rename('n_sampled_frames').reset_index()
        hbtime['occupancy_percent']=100*hbtime.n_sampled_frames/20
    else:
        hbocc=pd.DataFrame(columns=['hydrogen_bond','n_sampled_frames','occupancy_percent'])
        hbtime=pd.DataFrame(columns=['hydrogen_bond','ns_segment','n_sampled_frames','occupancy_percent'])
    hbocc.to_csv(odir/f'{label}_Hbond_occupancy.csv',index=False)
    hbtime.to_csv(odir/f'{label}_Hbond_1ns_timeline.csv',index=False)
    # Individual publication-ready graphs.
    fig,axs=plt.subplots(3,1,figsize=(11,10),sharex=True)
    _plot_curve(axs[0],df.time_ns,df.ligand_RMSD_vs_production_A,
                f'{label}: protein-aligned ligand heavy-atom displacement','RMSD vs first production frame (Å)')
    _plot_curve(axs[1],df.time_ns,df.ligand_RMSD_vs_original_docking_A,
                'Pose evolution relative to original docked coordinates','RMSD vs original docking pose (Å)')
    _plot_curve(axs[2],df.time_ns,df.Fe_ligand_min_heavy_A,
                'Closest heavy-atom approach to heme iron','Fe–ligand minimum distance (Å)')
    _save_png(fig,odir/f'{label}_pose_evolution.png')
    fig,axs=plt.subplots(2,1,figsize=(11,7),sharex=True)
    for col in ['ligand_bound_SASA_A2','ligand_isolated_SASA_A2','ligand_buried_SASA_A2']:
        axs[0].plot(ds.time_ns,ds[col],label=col.replace('_',' ').replace(' A2',''),lw=1)
    axs[0].legend(fontsize=8);axs[0].set_ylabel('Ligand SASA (Å²)');axs[0].set_title(f'{label}: ligand solvent exposure / burial')
    axs[1].plot(ds.time_ns,ds.ligand_bound_NOS_SASA_A2,label='N/O/S heavy atoms',lw=1)
    axs[1].plot(ds.time_ns,ds.ligand_bound_other_heavy_SASA_A2,label='Other heavy atoms',lw=1)
    axs[1].legend(fontsize=8);axs[1].set_ylabel('Bound ligand SASA (Å²)');axs[1].set_xlabel('Time (ns)')
    _save_png(fig,odir/f'{label}_SASA.png')
    fig,axs=plt.subplots(2,1,figsize=(11,7))
    axs[0].plot(hbc.time_ns,hbc.n_hbonds,ls='none',marker='.',markersize=1.5,alpha=0.6)
    axs[0].set(ylabel='Detected direct H-bonds',title=f'{label}: 50 ns ligand–protein H-bond dynamics')
    if not hbtime.empty:
        top_hb=hbocc.head(8).hydrogen_bond.tolist()
        matrix=np.zeros((len(top_hb),50),dtype=float)
        for k,name in enumerate(top_hb):
            subset=hbtime[hbtime.hydrogen_bond==name]
            for row in subset.itertuples(): matrix[k,int(row.ns_segment)-1]=row.occupancy_percent
        img=axs[1].imshow(matrix,aspect='auto',origin='upper',extent=[.5,50.5,len(top_hb)-.5,-.5],vmin=0,vmax=100,cmap='viridis')
        axs[1].set_yticks(np.arange(len(top_hb)));axs[1].set_yticklabels(top_hb,fontsize=7)
        fig.colorbar(img,ax=axs[1],label='Occupancy per 1 ns (%)')
    axs[1].set_xlabel('Production time (ns)');axs[1].set_title('Top hydrogen bonds by time interval')
    _save_png(fig,odir/f'{label}_Hbond_dynamics.png')
    allocc=contact_occup[contact_occup.contact_type=='all_heavy'].head(12)
    top_contacts=allocc.residue.tolist()
    ctime=contact_timeline[contact_timeline.contact_type=='all_heavy']
    mat=np.zeros((len(top_contacts),50),dtype=float)
    for k,name in enumerate(top_contacts):
        rows=ctime[ctime.residue==name]
        for r in rows.itertuples():mat[k,int(r.ns_segment)-1]=r.occupancy_percent
    fig,ax=plt.subplots(figsize=(11,6))
    im=ax.imshow(mat,aspect='auto',vmin=0,vmax=100,origin='upper',extent=[.5,50.5,len(top_contacts)-.5,-.5],cmap='viridis')
    ax.set_yticks(np.arange(len(top_contacts)));ax.set_yticklabels(top_contacts,fontsize=9)
    ax.set(xlabel='Production time (ns)',title=f'{label}: top residue contacts, 4.5 Å')
    fig.colorbar(im,ax=ax,label='Fraction of frames in each 1 ns interval (%)')
    _save_png(fig,odir/f'{label}_contact_time_heatmap.png')
    fig,axs=plt.subplots(2,1,figsize=(11,7),sharex=False)
    axs[0].plot(df.time_ns,df.Fe_S_A,lw=.7)
    axs[0].set(xlabel='Time (ns)',ylabel='Fe–S (Å)',title=f'{label}: heme–Cys coordination')
    axs[1].plot(rmsf.pdb_resid,rmsf.CA_RMSF_A,lw=0.8)
    axs[1].set(xlabel='Original PDB residue number',ylabel='Cα RMSF (Å)',title='Protein Cα mobility')
    _save_png(fig,odir/f'{label}_heme_and_RMSF.png')
    if docking['stereo_bonds']:
        fig,ax=plt.subplots(figsize=(10,3.5))
        for k,b in enumerate(docking['stereo_bonds']):
            ax.plot(df.time_ns,df[f'doublebond_{k+1}_dihedral_deg'],lw=.8,
                    label=f'3D SDF {b["reference_geometry"]}; original {b["original_dihedral_deg"]:.1f}°')
        ax.set(xlabel='Time (ns)',ylabel='Dihedral angle (degrees)',title=f'{label}: substituted alkene geometry')
        ax.set_ylim(-190,190);ax.legend(fontsize=8)
        _save_png(fig,odir/f'{label}_alkene_dihedral.png')
    # Export first, median-representative (10–50 ns), final structures for molecular viewers.
    rrow=df.loc[df.time_ns.between(10,50)].copy()
    med=float(rrow.ligand_RMSD_vs_original_docking_A.median())
    ridx=(rrow.ligand_RMSD_vs_original_docking_A-med).abs().idxmin()
    representative=df.loc[ridx]
    snapshot_rows={
        'first':df.iloc[0],
        'median_pose':representative,
        'final':df.iloc[-1],
    }
    exports={}
    for tag,row in snapshot_rows.items():
        ns=int(row.ns_segment);fi=int(row.frame_in_segment)
        exports[tag]=save_snapshot(label,tag,seg_lookup[ns][0],fi,client,bucket_name,
                                    top,info,docking['mol'],docking['atom_map'],reference,odir)
    meta['representative_frames']={tag:{'time_ns':float(row.time_ns),'ns_segment':int(row.ns_segment),
                                             'frame_in_segment':int(row.frame_in_segment),**exports[tag]}
                                   for tag,row in snapshot_rows.items()}
    # Ready-to-run PyMOL rendering recipe for a 3D binding-site illustration.
    pml = odir / f'{label}_publication_pose.pml'
    pml.write_text(
        f'load {exports["median_pose"]["pdb"]}, complex\n'
        'hide everything, all\n'
        'show cartoon, polymer.protein\n'
        'color gray70, polymer.protein\n'
        'show sticks, resn LIG or resn HEM\n'
        'color yellow, resn LIG\n'
        'color salmon, resn HEM\n'
        'select nearby, byres (polymer.protein within 4.5 of resn LIG)\n'
        'show sticks, nearby\n'
        'color cyan, nearby\n'
        'bg_color white\n'
        'orient resn LIG\n'
        'set ray_opaque_background, 0\n'
        f'png {label}_binding_site_3D.png, width=1800, height=1400, dpi=300, ray=1\n',
        encoding='utf-8')
    meta['pymol_render_script'] = pml.name
    meta['results']={
        'n_analyzed_frames':N,
        'n_sasa_frames':len(ds),
        'n_hbond_sampled_frames':total_hb_sampled,
        'ligand_RMSD_original_docking_mean_A':float(df.ligand_RMSD_vs_original_docking_A.mean()),
        'ligand_RMSD_firstproduction_mean_A':float(df.ligand_RMSD_vs_production_A.mean()),
        'bound_ligand_SASA_mean_A2':float(ds.ligand_bound_SASA_A2.mean()),
        'buried_ligand_SASA_mean_A2':float(ds.ligand_buried_SASA_A2.mean()),
        'direct_ligand_protein_Hbond_mean_per_sample':float(hbc.n_hbonds.mean()),
        'mean_Fe_S_A':float(df.Fe_S_A.mean()),
        'top_direct_hbonds':hbocc.head(8).to_dict(orient='records'),
        'top_contacts':allocc.head(12).to_dict(orient='records'),
    }
    if docking['stereo_bonds']:
        meta['alkene_geometry']={}
        for k,b in enumerate(docking['stereo_bonds']):
            angles=np.radians(df[f'doublebond_{k+1}_dihedral_deg'].to_numpy())
            orig=float(np.radians(b['original_dihedral_deg']))
            # Consistent relative geometry: compare cosine side of the original dihedral;
            # not a full CIP stereo reassignment on every trajectory frame.
            same=np.mean(np.sign(np.cos(angles))==np.sign(np.cos(orig)))
            meta['alkene_geometry'][f'alkene_{k+1}']={'SDF_assignment':b['reference_geometry'],
                 'original_dihedral_deg':b['original_dihedral_deg'],
                 'fraction_frames_same_geometric_side':float(same)}
    (odir/f'{label}_analysis_summary.json').write_text(json.dumps(meta,indent=2),encoding='utf-8')
    (odir/f'{label}_analysis_README.txt').write_text(
        f'{label}: 50 ns standalone trajectory analysis\n'
        'The analysis uses 5000 frames across 50 1-ns DCD segments; ligand pose RMSD is '
        'protein-CA aligned against (i) first production frame and (ii) original SDF/PDB docking pose.\n'
        'SASA uses heavy-atom Shrake-Rupley radii with 80 sphere points, 4 snapshots/ns. '
        'Ligand isolated SASA uses the SAME snapshot conformations, removing the receptor; '
        'ligand buried SASA = isolated minus bound. Explicit water is omitted in SASA surfaces.\n'
        'Contacts use 4.5 A ligand-protein heavy atom proximity over all 5000 frames. '
        'C/S packing proxy is a geometric heuristic, not an energy estimate.\n'
        'Wernet-Nilsson direct ligand-protein H-bonds use every fifth frame, excluding water.\n'
        'SDF stereo assignments are inferred from original three-dimensional reference. '
        'Alkene dihedral traces are from the corresponding topological atoms in the MD coordinates.\n'
        'Representative PDB/SDF files are protein-aligned simulation snapshots. '
        'Use a molecular viewer (PyMOL/ChimeraX) for publication 3D poses.\n'
        'All figures describe one trajectory per complex, not replicate variability.\n',encoding='utf-8')
    print(f'  {label}: advanced analysis finished. Dock RMSD {meta["results"]["ligand_RMSD_original_docking_mean_A"]:.2f} Å; '
          f'bound ligand SASA {meta["results"]["bound_ligand_SASA_mean_A2"]:.1f} Å²')
    return meta


def run_all():
    client,bucket_name=get_s3_client()
    selections=list_s3_objects(client,bucket_name)
    outcomes={}
    for label in LABELS:
        outcomes[label]=run_one(label,selections[label],client,bucket_name)
    (OUTPUT/'run_summary.json').write_text(json.dumps({label:outcomes[label]['results'] for label in LABELS},indent=2),encoding='utf-8')
    zpath=OUTPUT.parent/'CAL028_CAL015_advanced_analysis_results.zip'
    with zipfile.ZipFile(zpath,'w',compression=zipfile.ZIP_DEFLATED,compresslevel=6) as z:
        for f in sorted(OUTPUT.rglob('*')):
            if f.is_file() and f.name!='_temp_segment.dcd':
                z.write(f,arcname=str(f.relative_to(OUTPUT)))
    print('\nDONE. Download results:',zpath,'\nNo MD performed, original S3 files unchanged.')
    return zpath

if __name__=='__main__':
    run_all()
