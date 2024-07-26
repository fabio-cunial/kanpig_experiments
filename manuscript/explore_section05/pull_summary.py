import json
import pandas as pd
import joblib
import pysam
from collections import Counter
import ast

def add_summary(data):
    data['as_compP'] = data['as_TP'] + data['as_FP']
    data['as_baseP'] = data['as_TP'] + data['as_FN']
    data['as_compN'] = data['as_TN'] + data['as_FN']
    data['as_baseN'] = data['as_TN'] + data['as_FP']
    data['as_ppv'] = data['as_TP'] / data['as_compP']
    data['as_tpr'] = data['as_TP'] / data['as_baseP']
    data['as_tnr'] = data['as_TN'] / data['as_baseN']
    data['as_npv'] = data['as_TN'] / data['as_compN']
    data['as_acc'] = (data['as_TP'] + data['as_TN']) / (data['as_baseP'] + data['as_baseN'])
    data['as_ba'] = (data['as_tpr'] + data['as_tnr']) / 2
    data['as_f1'] = 2 * ((data['as_ppv'] * data['as_tpr']) / (data['as_ppv'] * data['as_tpr']))

def parse(name):
    with open(name + '/summary.json', 'r') as fh:
        d = json.load(fh)

    d['program'] = name.split('_')[1]

    d['gt_matrix']['(0, 0)'] = Counter()
    d['Concordant'] = d['TP-comp_TP-gt']
    d['Discordant'] = d['TP-comp_FP-gt']
    d['Missing'] = 0
    d['Present FP'] = 0
    d['Present TP'] = 0
    v = pysam.VariantFile(name + "/fp.vcf.gz")
    for entry in v:
        cgt = entry.samples[0]['GT']
        d['gt_matrix']['(0, 0)'][str(cgt)] += 1
        if None in cgt:
            d['Missing'] += 1
        elif 1 in cgt: # Only penalize for non-refhom
            d['Discordant'] += 1
            d['Present FP'] += 1
    # reclassify discordant which should be considered missing
    v = pysam.VariantFile(name + '/tp-comp.vcf.gz')
    for entry in v:
        cgt = entry.samples[0]['GT']
        if None in cgt:
            d['Missing'] += 1
            d['Discordant'] -= 1
        elif 1 in cgt:
            d['Present TP'] += 1
    # And by counting the Present TP/FP, we can re-estimate the precision/f1/recall of the variants
    d['gt_accuracy'] = d['gt_concordance']
    d['gt_concordance'] = d['Concordant'] / (d['Concordant'] + d['Discordant'])
    
    astate = []
    d['as_TP'] = 0
    d['as_TN'] = 0
    d['as_FP'] = 0
    d['as_FN'] = 0
    # Shouldn't happen for assembly vcfs
    d['as_filtered'] = 0
    d['as_missing'] = 0
    for b_key, b_data in d['gt_matrix'].items():
        b = list(ast.literal_eval(b_key))
        b = sorted([_ if _ is not None else -1 for _ in b])
        for c_key, c_data in b_data.items():
            c = list(ast.literal_eval(c_key))
            c = sorted([_ if _ is not None else -1 for _ in c])
            for i,j in zip(b, c):
                if i == -1:
                    d['as_filtered'] += c_data
                elif j == -1:
                    d['as_missing'] += c_data
                else:
                    state = "T" if i == j else "F"
                    condition = "P" if j == 1 else "N"
                    d[f'as_{state + condition}'] += c_data

    # Cleaning
    del(d['gt_matrix'])
    dkeys = [_ for _ in d.keys() if _.endswith('-gt')]
    for i in dkeys:
        del(d[i])
    return d

parts = ["bench_orig",
         "bench_kanpig",
         "bench_svjedi",
         "bench_snif",
         "bench_cutesv",
         ]
tmp = []
for i in parts:
    tmp.append(parse(i))

out = pd.DataFrame(tmp)
add_summary(out)
joblib.dump(out, 'results.jl')
print(out)
