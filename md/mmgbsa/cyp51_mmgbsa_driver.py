#!/usr/bin/env python3
"""CYP51 CAL028/CAL015 single-trajectory MM/GBSA, AmberTools / Google Drive.
Read-only access to MD source files; output written to a separate directory.

Invocation (in the AmberTools conda environment):
  python cyp51_mmgbsa_driver.py --system CAL028 --mode smoke
  python cyp51_mmgbsa_driver.py --system CAL028 --mode full
  python cyp51_mmgbsa_driver.py --system CAL015 --mode full
  python cyp51_mmgbsa_driver.py --system CAL028 --mode decomp
  python cyp51_mmgbsa_driver.py --system CAL028 --mode report

No MM/GBSA values are predicted by this script: all values are parsed from
successful AmberTools MMPBSA.py output on the supplied DCD trajectories.
"""
import argparse
import csv
import json
import os
import re
import shutil
import subprocess
import sys
from pathlib import Path
from datetime import datetime, timezone
import numpy as np
import parmed as pmd

DEFAULT_SOURCE = Path('/content/drive/MyDrive/CYP51_MD_animation_sources')
DEFAULT_OUTPUT = Path('/content/drive/MyDrive/CYP51_MD_MMGBSA')

SOLVENT_IONS = {
    'WAT','HOH','SOL','TIP3','TIP3P','TIP4','TIP4P','OPC','OPC3','SPC','SPCE',
    'NA','Na+','SOD','CL','Cl-','CLA','K','K+','POT','MG','MG2','Mg2+','CA','Ca2+',
}

TERMS = ['VDWAALS','EEL','EGB','ESURF','DELTA G gas','DELTA G solv','DELTA TOTAL']


def exec_cmd(args, log_path, cwd=None):
    args = [str(x) for x in args]
    log_path = Path(log_path)
    log_path.parent.mkdir(parents=True, exist_ok=True)
    print('$', ' '.join(args), flush=True)
    with log_path.open('w') as log:
        p = subprocess.Popen(args, cwd=cwd, stdout=subprocess.PIPE,
                             stderr=subprocess.STDOUT, text=True, bufsize=1)
        for line in p.stdout:
            log.write(line)
            log.flush()
            if ('WARNING' in line.upper() or 'ERROR' in line.upper() or
                'Processing' in line or 'Timing' in line or 'MMPBSA' in line):
                print(line.rstrip()[:250], flush=True)
        code = p.wait()
    if code:
        tail='\n'.join(log_path.read_text(errors='replace').splitlines()[-45:])
        raise RuntimeError(f'Command returned {code}. Log: {log_path}\n{tail}')
    return log_path


def audit_topology(source):
    parm = pmd.load_file(str(source))
    names = [r.name for r in parm.residues]
    ligand_res = [r for r in parm.residues if r.name == 'LIG']
    heme_res = [r for r in parm.residues if r.name == 'HEM']
    if len(ligand_res)!=1 or len(heme_res)!=1:
        raise RuntimeError(f'Expected exactly one LIG and one HEM residue; found '
                           f'{len(ligand_res)} LIG, {len(heme_res)} HEM. '
                           'Do not guess ligand masks.')
    cross_bonds=[]
    for b in parm.bonds:
        a1,a2=b.atom1,b.atom2
        if (a1.residue.name == 'LIG') != (a2.residue.name == 'LIG'):
            cross_bonds.append((str(a1), str(a2)))
    if cross_bonds:
        raise RuntimeError(f'Ligand has cross-component covalent bonds: {cross_bonds[:5]}. '
                           'Standard noncovalent MM/GBSA splitting is invalid.')
    heme_cys = [(str(b.atom1),str(b.atom2)) for b in parm.bonds
                if ((b.atom1.residue.name == 'HEM') != (b.atom2.residue.name == 'HEM'))
                and ('HEM' in (b.atom1.residue.name, b.atom2.residue.name))]
    strip_names = sorted(set(names) & SOLVENT_IONS)
    if not strip_names:
        raise RuntimeError('No water/ion residues identified; inspect the topology before stripping.')
    expected_solute_atoms = sum(len(r.atoms) for r in parm.residues if r.name not in strip_names)
    # Radius provenance: do not silently assume correct Fe solvation physics.
    fe = [a for a in parm.atoms if a.residue.name=='HEM' and a.name.upper()=='FE']
    fe_radii = [getattr(a, 'solvent_radius', None) for a in fe]
    meta = {
        'original_atoms': len(parm.atoms),
        'original_residues': len(parm.residues),
        'ligand_resid': ligand_res[0].idx + 1,
        'ligand_atom_count': len(ligand_res[0].atoms),
        'ligand_net_charge_e': float(sum(a.charge for a in ligand_res[0].atoms)),
        'heme_resid': heme_res[0].idx + 1,
        'heme_protein_cross_bonds': heme_cys,
        'heme_FE_solvation_radii_A': fe_radii,
        'strip_residue_names': strip_names,
        'expected_dry_atoms': expected_solute_atoms,
        'output_format': 'kcal/mol',
        'notes': 'Heme included in the receptor; ligand remains its original charged form.'
    }
    print('Source topology:', source)
    print('Atoms:', len(parm.atoms), 'Dry complex atoms:', expected_solute_atoms)
    print('LIG residue:', meta['ligand_resid'], 'charge:', round(meta['ligand_net_charge_e'], 4))
    print('HEM residue:', meta['heme_resid'], 'Fe radii:', fe_radii)
    print('Heme-protein bonds:', heme_cys[:4])
    print('Water/ion labels to strip:', strip_names)
    return meta


def prepare_topologies(system, source_root, out_root):
    source = source_root / system / 'system.prmtop'
    if not source.is_file():
        raise FileNotFoundError(source)
    out = out_root/system
    out.mkdir(parents=True, exist_ok=True)
    auditfile=out/'topology_audit.json'
    meta = audit_topology(source)
    auditfile.write_text(json.dumps(meta,indent=2))
    complexp,receptorp,ligandp=[out/name for name in ('complex.prmtop','receptor.prmtop','ligand.prmtop')]
    if not all(p.is_file() and p.stat().st_size>1000 for p in (complexp,receptorp,ligandp)):
        mask=':'+','.join(meta['strip_residue_names'])
        cmd=['ante-MMPBSA.py', '-p', source, '-c', complexp, '-r', receptorp,
             '-l', ligandp, '-s', mask, '-n', ':LIG', '--radii=mbondi2']
        exec_cmd(cmd,out/'ante_mmpbsa.log',cwd=out)
    comp,rec,lig=(pmd.load_file(str(p)) for p in (complexp,receptorp,ligandp))
    if len(comp.atoms)!=meta['expected_dry_atoms']:
        raise RuntimeError(f'Dry topology atom count {len(comp.atoms)} does not match '
                           f'expected {meta["expected_dry_atoms"]}. Check solvent mask.')
    if len(rec.atoms)+len(lig.atoms)!=len(comp.atoms):
        raise RuntimeError('Receptor + ligand atom totals differ from complex.')
    if any(r.name=='HEM' for r in lig.residues) or not any(r.name=='HEM' for r in rec.residues):
        raise RuntimeError('HEM must remain in receptor, not the ligand.')
    if abs(sum(a.charge for a in lig.atoms)-meta['ligand_net_charge_e'])>0.001:
        raise RuntimeError('Ligand charge changed while preparing MM/GBSA topology.')
    # Some heme force fields produce missing PB radii, which should not be silently fixed.
    bad=[]
    for a in comp.atoms:
        radius=getattr(a,'solvent_radius',None)
        if radius is None or not np.isfinite(float(radius)) or float(radius)<=0:
            bad.append((a.residue.name, a.name, radius))
    if bad:
        raise RuntimeError(f'Invalid GB/SA radii for {len(bad)} atoms (first 12: '
                           f'{bad[:12]}). Need validate heme Fe/nonstandard atom radii '
                           'rather than claim an energy calculation with invalid parameters.')
    (out/'topology_audit.json').write_text(json.dumps({**meta,
        'dry_complex_atoms':len(comp.atoms), 'receptor_atoms':len(rec.atoms),
        'ligand_atoms':len(lig.atoms),'ligand_charge_after_strip':float(sum(a.charge for a in lig.atoms)),
        'GB_model':'igb=5 (OBC2), mbondi2 radii, saltcon=0.150 M'},indent=2))
    print(f'{system}: topology preparation and component checks PASSED.',flush=True)
    return meta


def block_dry(system, block, meta, source_root, out_root, smoke=False):
    out=out_root/system/('smoke' if smoke else f'block_{block:02d}')
    out.mkdir(parents=True,exist_ok=True)
    dry=out/'sample_dry.nc'
    if dry.exists() and dry.stat().st_size>10000:
        print(system, 'reusing prepared dry trajectory:', dry)
        return dry
    src=source_root/system
    begin=1 if smoke else (block-1)*10+1
    end=2 if smoke else block*10
    step=50 if smoke else 25   # 2 frames/ns test; 4 frames/ns for 50ns analysis
    lines=[f'parm {src/"system.prmtop"}']
    for ns in range(begin,end+1):
        dcd=src/f'ns_{ns:03d}.dcd'
        if not dcd.is_file(): raise FileNotFoundError(dcd)
        lines.append(f'trajin {dcd} 1 last {step}')
    lines += ['autoimage', 'strip :'+','.join(meta['strip_residue_names']),
              f'trajout {dry} netcdf', 'run', 'quit']
    script=out/'strip_and_sample.cpptraj'
    script.write_text('\n'.join(lines)+'\n')
    exec_cmd(['cpptraj','-i',script],out/'cpptraj.log',cwd=out)
    if not dry.exists() or dry.stat().st_size<10000:
        raise RuntimeError(f'cpptraj did not produce a valid-looking NetCDF: {dry}')
    # independent topology/trajectory validation with cpptraj itself
    check=out/'audit_dry.cpptraj'
    check.write_text(f'parm {out_root/system/"complex.prmtop"}\n'
                     f'trajin {dry} 1 1\nrun\nquit\n')
    exec_cmd(['cpptraj','-i',check],out/'audit_dry.log',cwd=out)
    print(system, 'prepared dry trajectory:',dry, flush=True)
    return dry


def gb_input(decomp=False):
    general='''&general\n startframe=1,\n endframe=999999,\n interval=%s,\n keep_files=0,\n verbose=1,\n/\n''' % (4 if decomp else 1)
    gb='''&gb\n igb=5,\n saltcon=0.150,\n/\n'''
    decompblock='''&decomp\n idecomp=2,\n dec_verbose=1,\n print_res="within 5",\n/\n''' if decomp else ''
    return general+gb+decompblock


def parse_results(p):
    s=Path(p).read_text(errors='replace')
    if 'INCONSISTENCIES EXIST' in s.upper():
        raise RuntimeError(f'MMPBSA reported internal-energy inconsistencies. Inspect {p}')
    if 'GENERALIZED BORN:' not in s:
        raise RuntimeError(f'No GENERALIZED BORN section found in {p}.')
    gb=s.split('GENERALIZED BORN:',1)[1].split('POISSON BOLTZMANN:',1)[0]
    anchor='Differences (Complex - Receptor - Ligand):'
    if anchor not in gb:
        raise RuntimeError(f'No DELTA component section found in {p}.')
    delta=gb.split(anchor,1)[1]
    result={}
    num=r'[-+]?(?:\d*\.\d+|\d+\.?\d*)(?:[Ee][-+]?\d+)?'
    for name in TERMS:
        m=re.search(r'^\s*'+re.escape(name)+r'\s+('+num+r')\s+('+num+r')',delta,re.M)
        if m:
            result[name]={'mean':float(m.group(1)),'sd':float(m.group(2))}
    if 'DELTA TOTAL' not in result:
        raise RuntimeError(f'Missing delta total in {p}. The source file is retained.')
    return result


def run_one(system, block, meta, source_root, out_root, smoke=False, decomp=False):
    bdir=out_root/system/('smoke' if smoke else f'block_{block:02d}')
    target_name = 'MMGBSA_DECOMP' if decomp else 'MMGBSA'
    target=bdir/f'{target_name}_summary.json'
    if target.is_file() and not smoke:
        try:
            data=json.loads(target.read_text())
            if 'DELTA TOTAL' in data['terms']:
                print(f'{system} block {block}: already completed, skipping.', flush=True)
                return data
        except (ValueError,KeyError): pass
    dry=block_dry(system,block,meta,source_root,out_root,smoke=smoke)
    infile=bdir/('mmpbsa_decomp.in' if decomp else 'mmpbsa.in')
    infile.write_text(gb_input(decomp=decomp))
    complexp,receptorp,ligandp=[out_root/system/name for name in ('complex.prmtop','receptor.prmtop','ligand.prmtop')]
    prefix='DECOMP' if decomp else 'GBSA'
    outfile=bdir/f'{prefix}_RESULTS.dat'
    cmd=['MMPBSA.py','-O','-i',infile,'-cp',complexp,'-rp',receptorp,'-lp',ligandp,
         '-y',dry, '-o',outfile, '-eo',bdir/f'{prefix}_perframe.csv']
    if decomp:
        cmd+=['-do',bdir/'DECOMP_RESIDUES.dat','-deo',bdir/'DECOMP_RESIDUES.csv']
    print(f'\n=== {system}: {"SMOKE" if smoke else f"{(block-1)*10}-{block*10} ns"}'
          f' (GB-OBC2, {"decomposition" if decomp else "standard"}) ===',flush=True)
    scratch=Path('/content/CYP51_MMGBSA_scratch')/system/('smoke' if smoke else f'block_{block:02d}')
    scratch.mkdir(parents=True,exist_ok=True)
    exec_cmd(cmd,bdir/f'{prefix}_run.log',cwd=scratch)
    terms=parse_results(outfile)
    result={'system':system,'block':block,'smoke':smoke,'decomposition':decomp,
            'frames_expected':(4 if smoke else 10 if decomp else 40),
            'time_ns':[0,2] if smoke else [(block-1)*10,block*10],
            'terms':terms,'result_file':str(outfile),'trajectory':str(dry),
            'calculation_model':'GB-OBC2 (igb 5), mbondi2, saltcon 0.15 M; no entropy term',
            'datetime_utc':datetime.now(timezone.utc).isoformat()}
    target.write_text(json.dumps(result,indent=2))
    print(f'{system} block {block} ΔG_MMGBSA (no entropy): '
          f'{terms["DELTA TOTAL"]["mean"]:.3f} ± {terms["DELTA TOTAL"]["sd"]:.3f} kcal/mol',flush=True)
    return result


def summarize(system,out_root):
    """Create block summaries, including partial results; plotting is optional."""
    import csv
    from statistics import stdev
    folder=out_root/system
    results=[]
    for block in range(1,6):
        p=folder/f'block_{block:02d}'/'MMGBSA_summary.json'
        if p.is_file():
            results.append(json.loads(p.read_text()))
    if not results:
        print(f'{system}: no completed 10ns blocks yet.')
        return
    records=[]
    for r in sorted(results,key=lambda x:x['block']):
        record={'system':system,'block':r['block'],
                'midpoint_ns':(r['block']-0.5)*10,
                'n_frames':r['frames_expected']}
        for key,value in r['terms'].items():
            record[key]=value['mean']
        records.append(record)
    fields=['system','block','midpoint_ns','n_frames'] + TERMS
    with (folder/'MMGBSA_10ns_block_means.csv').open('w',newline='') as handle:
        writer=csv.DictWriter(handle,fieldnames=fields,extrasaction='ignore')
        writer.writeheader();writer.writerows(records)
    frames=sum(r['n_frames'] for r in records)
    pooled={key:(sum(r[key]*r['n_frames'] for r in records if key in r) /
                 sum(r['n_frames'] for r in records if key in r))
            for key in TERMS if any(key in r for r in records)}
    means=[r['DELTA TOTAL'] for r in records]
    pooled.update({'n_completed_10ns_blocks':len(records),
                   'n_analyzed_frames':frames,
                   'block_to_block_SD_kcal_mol':stdev(means) if len(means)>1 else None,
                   'block_mean_SEM_kcal_mol':stdev(means)/len(means)**0.5 if len(means)>1 else None,
                   'full_50ns_covered':len(records)==5})
    (folder/'MMGBSA_aggregate_summary.json').write_text(json.dumps(pooled,indent=2))
    try:
        import matplotlib
        matplotlib.use('Agg')
        import matplotlib.pyplot as plt
        fig,ax=plt.subplots(figsize=(8,4.7))
        ax.plot([r['midpoint_ns'] for r in records],means,'o-',linewidth=1.8)
        ax.set(xlim=(0,50),xlabel='Simulation time (ns)',
               ylabel=r'MM/GBSA $\Delta G$ (kcal/mol)',
               title=f'{system}: MM/GBSA by 10 ns block')
        ax.grid(alpha=0.2)
        fig.tight_layout()
        fig.savefig(folder/'MMGBSA_10ns_block_trend.png',dpi=300)
        plt.close(fig)
    except ImportError:
        print('Plot skipped (matplotlib unavailable); numerical reports are saved.',flush=True)
    print(f'{system}: {len(records)} blocks, {frames} frames reported.',flush=True)
    if len(records)<5:
        print('PARTIAL RESULT — 50 ns analysis not complete.',flush=True)

def main():
    ap=argparse.ArgumentParser()
    ap.add_argument('--system',required=True,choices=['CAL028','CAL015'])
    ap.add_argument('--mode',required=True,choices=['smoke','full','report','decomp'])
    ap.add_argument('--source-root',default=str(DEFAULT_SOURCE))
    ap.add_argument('--out-root',default=str(DEFAULT_OUTPUT))
    args=ap.parse_args()
    source_root=Path(args.source_root)
    out_root=Path(args.out_root)
    if args.mode=='report':
        summarize(args.system,out_root);return
    for binary in ('ante-MMPBSA.py','cpptraj','MMPBSA.py'):
        if not shutil.which(binary):
            raise RuntimeError(f'{binary} not found in PATH. Run in the dedicated AmberTools environment.')
    meta=prepare_topologies(args.system,source_root,out_root)
    if args.mode=='smoke':
        run_one(args.system,1,meta,source_root,out_root,smoke=True)
    elif args.mode=='full':
        for block in range(1,6):
            run_one(args.system,block,meta,source_root,out_root)
        print(f'{args.system}: all five blocks completed; run --mode report separately.', flush=True)
    elif args.mode=='decomp':
        for block in range(1,6):
            run_one(args.system,block,meta,source_root,out_root,decomp=True)
    print('Done. Source DCDs/topologies were not modified.',flush=True)

if __name__=='__main__':
    main()
