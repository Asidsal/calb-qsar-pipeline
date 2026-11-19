"""Offline report recovery for completed CAL028/CAL015 MM/GBSA.

Reads only the previously computed 5 x 10 ns block outputs for each compound.
Does NOT run MD, cpptraj, AmberTools, or MMPBSA.py.

Run in Google Colab: %run /content/CYP51_MMGBSA_Report_Recovery.py
Or copy this file into Colab and invoke main().
"""
from pathlib import Path
import csv
import json
import math
import zipfile

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np

ROOT = Path('/content/drive/MyDrive/CYP51_MD_MMGBSA')
SYSTEMS = ('CAL028', 'CAL015')
COMPONENTS = ('VDWAALS', 'EEL', 'EGB', 'ESURF', 'DELTA G gas', 'DELTA G solv', 'DELTA TOTAL')


def write_csv(destination, rows, columns):
    with Path(destination).open('w', newline='', encoding='utf-8') as handle:
        writer = csv.DictWriter(handle, fieldnames=columns)
        writer.writeheader()
        writer.writerows(rows)


def calculate(root=ROOT):
    root = Path(root)
    if not root.is_dir():
        raise FileNotFoundError(f'Google Drive MMGBSA results folder missing: {root}')

    all_summaries = {}
    include = []
    for system in SYSTEMS:
        folder = root / system
        block_rows, data = [], {}
        for block in range(1, 6):
            block_folder = folder / f'block_{block:02d}'
            summary_path = block_folder / 'MMGBSA_summary.json'
            if not summary_path.is_file():
                raise FileNotFoundError(f'Missing completed calculation: {summary_path}')
            item = json.loads(summary_path.read_text(encoding='utf-8'))
            terms = item.get('terms', {})
            expected = int(item.get('frames_expected', 0))
            if int(item.get('block', block)) != block or expected != 40:
                raise ValueError(f'{summary_path}: unexpected block number or frame count')
            missing = [term for term in COMPONENTS if term not in terms]
            if missing:
                raise ValueError(f'{summary_path}: missing terms {missing}')
            data[block] = item
            row = {'system': system, 'block': block, 'start_ns': (block-1)*10,
                   'end_ns': block*10, 'midpoint_ns': (block-0.5)*10, 'n_frames': expected}
            for term in COMPONENTS:
                row[f'{term} mean'] = float(terms[term]['mean'])
                row[f'{term} sd'] = float(terms[term]['sd'])
            block_rows.append(row)
            include.append(summary_path)
            for basename in ('GBSA_RESULTS.dat', 'GBSA_perframe.csv', 'GBSA_run.log'):
                p = block_folder / basename
                if p.is_file():
                    include.append(p)

        means, sds, total_n, block_sd = {}, {}, 0, {}
        for term in COMPONENTS:
            vals = np.array([row[f'{term} mean'] for row in block_rows], dtype=float)
            deviations = np.array([row[f'{term} sd'] for row in block_rows], dtype=float)
            ns = np.array([row['n_frames'] for row in block_rows], dtype=int)
            N = int(ns.sum())
            mu = float(np.average(vals, weights=ns))
            # Reconstruct sample SD across all frames from group sample SDs and means.
            pooled_variance = (np.sum((ns-1)*deviations**2) + np.sum(ns*(vals-mu)**2)) / (N-1)
            means[term] = mu
            sds[term] = float(np.sqrt(max(pooled_variance, 0)))
            block_sd[term] = float(np.std(vals, ddof=1))
            total_n = N
        if total_n != 200:
            raise ValueError(f'{system}: expected exactly 200 evaluated frames, got {total_n}')

        aggregation = {
            'system': system,
            'sampling': '5 x 10 ns, 40 frames per block; 200 frames from 50 ns',
            'method': 'AmberTools MM/GBSA GB-OBC2 igb=5, mbondi2, saltcon=0.15 M, no entropy term',
            'n_analyzed_frames': total_n,
            'full_50ns_covered': True,
            'aggregate_terms': {term: {
                'mean': means[term], 'pooled_frame_sd': sds[term],
                'block_mean_sd': block_sd[term]}
                for term in COMPONENTS},
            'note': 'Pooled SD reconstructed from block mean and sample SD; rounded input may cause small differences from frame-level exact SD.'
        }

        out_json = folder / 'MMGBSA_aggregate_summary.json'
        out_json.write_text(json.dumps(aggregation, indent=2), encoding='utf-8')
        out_csv = folder / 'MMGBSA_10ns_block_means.csv'
        write_csv(out_csv, block_rows, list(block_rows[0].keys()))
        include.extend([out_json, out_csv])

        mid = [x['midpoint_ns'] for x in block_rows]
        y = [x['DELTA TOTAL mean'] for x in block_rows]
        err = [x['DELTA TOTAL sd'] for x in block_rows]
        fig, ax = plt.subplots(figsize=(7.1,4.3), dpi=160)
        ax.errorbar(mid, y, yerr=err, fmt='o-', color='black', capsize=3,
                    linewidth=1.4, markersize=4, ecolor='0.55')
        ax.axhline(means['DELTA TOTAL'], color='0.50', linestyle='--', linewidth=1.0,
                   label=f'50 ns mean: {means["DELTA TOTAL"]:.2f} kcal/mol')
        ax.set_xlabel('Simulation time (ns)')
        ax.set_ylabel('MM/GBSA ΔTOTAL (kcal/mol)')
        ax.set_title(f'{system}: MM/GBSA across 50 ns')
        ax.set_xlim(0,50)
        ax.set_xticks([0,10,20,30,40,50])
        ax.legend(frameon=False, loc='best', fontsize=8)
        ax.spines[['top','right']].set_visible(False)
        fig.tight_layout()
        plot = folder/'MMGBSA_10ns_block_trend.png'
        fig.savefig(plot, dpi=300, bbox_inches='tight')
        plt.close(fig)
        include.append(plot)

        terms_plot = ['VDWAALS','EEL','EGB','ESURF']
        fig, ax = plt.subplots(figsize=(7.3,4.3), dpi=160)
        xs = np.arange(len(terms_plot))
        components = np.array([means[t] for t in terms_plot])
        ax.bar(xs, components, color='0.45', edgecolor='black', linewidth=.75)
        ax.axhline(0, color='black', linewidth=.75)
        ax.set_xticks(xs, terms_plot)
        ax.set_ylabel('MM/GBSA energy (kcal/mol)')
        ax.set_title(f'{system}: MM/GBSA energy components')
        ax.spines[['top','right']].set_visible(False)
        fig.tight_layout()
        component_plot=folder/'MMGBSA_energy_components.png'
        fig.savefig(component_plot, dpi=300, bbox_inches='tight')
        plt.close(fig)
        include.append(component_plot)
        all_summaries[system]=aggregation
        print(f'{system}: VERIFIED {total_n}/200 frames; ΔTOTAL = '
              f'{means["DELTA TOTAL"]:.2f} ± {sds["DELTA TOTAL"]:.2f} kcal/mol')

    overview=[]
    for system, item in all_summaries.items():
        row={'system':system,'n_frames':item['n_analyzed_frames']}
        for term in COMPONENTS:
            row[f'{term} mean']=item['aggregate_terms'][term]['mean']
            row[f'{term} sd']=item['aggregate_terms'][term]['pooled_frame_sd']
        overview.append(row)
    overview_csv=root/'MMGBSA_both_compounds_summary.csv'
    write_csv(overview_csv, overview, list(overview[0].keys()))
    include.append(overview_csv)

    archive=root/'CYP51_MMGBSA_results_share.zip'
    with zipfile.ZipFile(archive,'w',compression=zipfile.ZIP_DEFLATED) as z:
        for p in include:
            if p.exists():
                z.write(p, p.relative_to(root))
    print(f'REPORT SUCCESSFUL: {archive} ({archive.stat().st_size/1048576:.2f} MiB)')
    print('Original MMPBSA results and trajectory data were not modified.')
    return archive


if __name__ == '__main__':
    calculate()
