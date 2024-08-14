"""
Given stats results from analysis.py, turn into giant flat tables

stats should have name SAMPLE.COVERAGE.TECH.stats.jl
Where SAMPLE and COVERAGE will be added as columns to the tables

Creates 4 tables: gt_dist.table.txt, svtype.table.txt, neigh.table.txt, intersect.table.txt
"""
import os
import sys
import joblib
import pandas as pd

# How to consolidate the gt distribution for a single joblib file

def make_frame(rows, sample, coverage, exp, tech, destination):
    if not rows:
        return
    frame = pd.concat(rows)
    frame['sample'] = sample
    frame['technology'] = tech
    frame['coverage'] = coverage
    frame['experiment'] = exp
    frame.reset_index(drop=True, inplace=True)
    destination.append(frame)

for subset in ['all', 'tr', 'non_tr']:
    gt_dist = []
    svtype = []
    neigh = []
    intersect = []

    for i in sys.argv[1:]:
        name = os.path.basename(i)
        sample, coverage, tech, _ = name.split('.')
        data = joblib.load(i)
        # This now needs to loop through trs and non-trs
        data = data[subset]
        # Extract baseline genotype dist
        base = data['base_gt_dist']
        base['program'] = 'baseline'
        base['sample'] = sample
        base['technology'] = tech
        base['coverage'] = coverage
        base['experiment'] = 'truth'
        gt_dist.append(base)
        
        for exp in ['ts', 'ds', 'tm', 'dm']:
            gt_rows = []
            ty_rows = []
            ne_rows = []
            in_rows = []
            for prog in data[exp]:
                d = data[exp][prog]['gt_dist']
                d['program'] = prog
                gt_rows.append(d)

                d = data[exp][prog]['type']
                d['program'] = prog
                ty_rows.append(d)

                d = data[exp][prog]['neigh']
                d['program'] = prog
                ne_rows.append(d)

                if exp != 'ts':
                    d = data[exp][prog]['table']
                    d['program'] = prog
                    in_rows.append(d)

            make_frame(gt_rows, sample, coverage, exp, tech, gt_dist)
            make_frame(ty_rows, sample, coverage, exp, tech, svtype)
            make_frame(ne_rows, sample, coverage, exp, tech, neigh)
            if exp != 'ts':
                make_frame(in_rows, sample, coverage, exp, tech, intersect)

    gt_dist = pd.concat(gt_dist)
    gt_dist.to_csv(f"{subset}.gt_dist.table.txt", sep='\t', index=False)

    svtype = pd.concat(svtype)
    svtype.to_csv(f"{subset}.svtype.table.txt", sep='\t', index=False)

    neigh = pd.concat(neigh)
    neigh.to_csv(f"{subset}.neighbor.table.txt", sep='\t', index=False)

    intersect = pd.concat(intersect)
    intersect.to_csv(f"{subset}.intersect.table.txt", sep='\t', index=False)
