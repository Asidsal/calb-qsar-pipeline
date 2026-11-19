"""Static validation of supplemental CYP51 MD data; no AmberTools required."""
from pathlib import Path
import ast
import csv
import json

ROOT = Path(__file__).resolve().parent.parent

def test_extension():
    for cmpd in ('CAL028','CAL015'):
        base = ROOT/'md/analysis_50ns/advanced'/cmpd
        assert base.is_dir()
        for kind, n in [('pose_dynamics_timeseries',5000),('SASA_timeseries',200),('Hbond_count_timeseries',1000)]:
            with (base/f'{cmpd}_{kind}.csv').open() as f:
                assert sum(1 for _ in csv.DictReader(f)) == n
        for kind in ['first_protein_heme_ligand','median_pose_protein_heme_ligand','final_protein_heme_ligand']:
            assert (base/f'{cmpd}_{kind}.pdb').is_file()
    s=(ROOT/'md/mmgbsa/cyp51_mmgbsa_driver.py').read_text()
    ast.parse(s)
    assert 'strip_mdcrd=0' not in s
    assert "summarize(args.system,out_root)" in s # invoked only by explicit report mode
    with (ROOT/'md/mmgbsa/results/MMGBSA_both_compounds_summary.csv').open() as f:
        records = {r['system']: r for r in csv.DictReader(f)}
    assert set(records) == {'CAL028','CAL015'}
    assert all(int(r['n_frames']) == 200 for r in records.values())
    assert abs(float(records['CAL028']['DELTA TOTAL mean']) + 60.42092) < 1e-5
    assert abs(float(records['CAL015']['DELTA TOTAL mean']) + 59.46166) < 1e-5
    for nm in ('CAL028','CAL015'):
        assert (ROOT/'md/visualization_3d'/f'{nm}_annotated_interactions.png').exists()
    assert not (ROOT/'manuscript').exists()  # publication files are maintained separately
    assert (ROOT/'docs/retrosynthesis/CAL028_proposed_convergent_synthesis.png').is_file()
    return True

if __name__=='__main__':
    assert test_extension()
    print('PASS: MD extension structure and source validation')
