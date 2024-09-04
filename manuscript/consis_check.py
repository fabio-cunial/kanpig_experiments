import sys
import pysam
import truvari
import numpy as np
from collections import Counter

header =["total", "pac_r9", "pac_r10", "r9_r10", "3x", "cor_pac_r9", "cor_pac_r10", "cor_r9_r10", "cor_3x"]
same_total = [0, 0, 0, 0, 0, 0, 0, 0, 0]
same_missing = [0, 0, 0, 0, 0, 0, 0, 0, 0]
same_present = [0, 0, 0, 0, 0, 0, 0, 0, 0]
SIZEMIN = 50
SIZEMAX = 10000

vcf = pysam.VariantFile(sys.argv[1])
bed = truvari.build_region_tree(vcfA=vcf, includebed=sys.argv[2])
vcf_i = truvari.region_filter(vcf, bed)

def update_counts(arr, gts):
    arr[0] += 1
    # pac9
    if gts[1] == gts[2]:
        arr[1] += 1
        if gts[0] == gts[1]:
            arr[5] += 1
    #pac10
    if gts[1] == gts[3]:
        if gts[0] == gts[1]:
            arr[6] += 1
        arr[2] += 1
    # 910
    if gts[2] == gts[3]:
        if gts[0] == gts[2]:
            arr[7] += 1
        arr[3] += 1
        #3x
        if gts[3] == gts[1]:
            arr[4] += 1
            if gts[0] == gts[1]:
                arr[8] += 1

for entry in vcf_i:
    sz = truvari.entry_size(entry)
    if sz < SIZEMIN or sz > SIZEMAX:
        continue
    if entry.chrom in ['chrX', 'chrY']:
        continue
    if None in entry.samples[0]['GT']:
        continue
    # just to see, must have 10x+ on all sites
    all_10x = True
    for samp in [1,2,3]:
        try:
            if entry.samples[samp]['SQ'] < 50:
                all_10x = False
        except Exception:
            all_10x = False
    if not all_10x:
        continue
    m_gts = []
    for samp in entry.samples:
        m_gts.append(truvari.get_gt(entry.samples[samp]['GT']))
    
    missing = truvari.GT.NON in m_gts
    any_ref = truvari.GT.REF in m_gts

    update_counts(same_total, m_gts)
    if not missing:
        update_counts(same_missing, m_gts)
    if not any_ref and not missing:
        update_counts(same_present, m_gts)

import pandas as pd

df = pd.DataFrame([same_total, same_missing, same_present], columns=header, index=['total', 'missing', 'present'])
print(df.to_csv(sep='\t'))
"""

def decode(key):
    val = []
    for i,v in enumerate(['base', 'pac', 'r9', 'r10']):
        if key & (1<<i):
            val.append(v)
        else:
            val.append(".")
    return "\t".join(val)
            
import json
print('total', total)
same = dict(same)
for key in sorted(same.keys()):
    v = same[key]
    print(decode(key), v, sep='\t')
# Just the three pairs:
import itertools
print("pairs")
for i,j in itertools.combinations((1,2,3), 2):
    tot = sum([same[k] for k in same if (k & (1<<i)) and (k & (1<<j))])
    print(i,j, tot, tot/total)

print("correct pairs")
for i,j in itertools.combinations((1,2,3), 2):
    tot = sum([same[k] for k in same if (k & 1) and (k & (1<<i)) and (k & (1<<j))])
    print(i,j, tot, tot/total)

print("same three")
tot = sum([same[k] for k in same if k in [14, 15]])
print(tot, tot / total)

print("same three correct")
tot = sum([same[k] for k in same if k == 15])
print(tot, tot / total)
"""



