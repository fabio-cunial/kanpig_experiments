"""
Given stats results from analysis.py, turn into giant flat tables

stats should have name SAMPLE.COVERAGE.stats.jl
Where SAMPLE and COVERAGE will be added as columns to the tables

Creates 4 tables: gt_dist.table.txt, svtype.table.txt, neigh.table.txt, intersect.table.txt
"""
import os
import sys
import joblib
import pandas as pd

# How to consolidate the gt distribution for a single joblib file
gt_dist = []
svtype = []
neigh = []
intersect = []

def make_frame(rows, sample, coverage, exp):
    frame = pd.concat(rows)
    frame['sample'] = sample
    frame['coverage'] = coverage
    frame['experiment'] = exp
    frame.reset_index(drop=True, inplace=True)
    return frame

for i in sys.argv[1:]:
    name = os.path.basename(i)
    sample, coverage, _, _ = name.split('.')
    data = joblib.load(i)

    # Extract baseline genotype dist
    base = data['base_gt_dist']
    base['program'] = 'baseline'
    base['coverage'] = coverage
    base['sample'] = sample
    base['experiment'] = 'truth'
    gt_dist.append(base)

    for exp in ['ts', 'ds', 'tm', 'dm']:
        gt_rows = []
        ty_rows = []
        ne_rows = []
        in_rows = []
        for key in data[exp]:
            d = data[exp][key]['gt_dist']
            d['program'] = key
            gt_rows.append(d)

            d = data[exp][key]['type']
            d['program'] = key
            ty_rows.append(d)

            d = data[exp][key]['neigh']
            d['program'] = key
            ne_rows.append(d)
            if exp != 'ts':
                d = data[exp][key]['table']
                d['program'] = key
                in_rows.append(d)

        gt_dist.append(make_frame(gt_rows, sample, coverage, exp))
        svtype.append(make_frame(ty_rows, sample, coverage, exp))
        neigh.append(make_frame(ne_rows, sample, coverage, exp))
        if exp != 'ts':
            intersect.append(make_frame(in_rows, sample, coverage, exp))

gt_dist = pd.concat(gt_dist)
gt_dist.to_csv("gt_dist.table.txt", sep='\t', index=False)

svtype = pd.concat(svtype)
svtype.to_csv("svtype.table.txt", sep='\t', index=False)

neigh = pd.concat(neigh)
neigh.to_csv("neighbor.table.txt", sep='\t', index=False)

intersect = pd.concat(intersect)
intersect.to_csv("intersect.table.txt", sep='\t', index=False)
