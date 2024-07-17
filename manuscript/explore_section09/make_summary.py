import sys
import pysam
import truvari
import json
import pandas as pd
import joblib

SIZEMIN = 50
SIZEMAX=10000
NEIGHLIMIT = 6
program='test'
technology='test'
# This should have numneigh annotation inside
base_vcf = sys.argv[1]
bed_fn = sys.argv[2]
comp_vcf = sys.argv[3]
bench_jl = sys.argv[4]

vcf = pysam.VariantFile(base_vcf)
bed = truvari.build_region_tree(vcfA=vcf, includebed=bed_fn)
if 'chrX' in bed:
    del(bed['chrX'])
if 'chrY' in bed:
    del(bed['chrY'])
for i in list(bed.keys()):
    if '_' in i:
        del(bed[i])

vcf_i = truvari.region_filter(vcf, bed)
m_cnt_base = {'Program':program,
         'Technology':technology,
         'Concordant': 0,
         'Discordant': 0,
         'Missing': 0,
         'Filtered': 0,
         'TP': 0, 'FP': 0, 'TN': 0, 'FN':0}
m_cnt = {'DEL': dict(m_cnt_base),
         'INS': dict(m_cnt_base),
         'TOT': dict(m_cnt_base)}

m_neigh_cnt = {}
for i in range(0, NEIGHLIMIT):
    m_neigh_cnt[i] = dict(m_cnt_base)

###

# Load the data
data = joblib.load(bench_jl)

# Separate out the variants from the base VCF and add new columns of the base/comp ids
base = data[data['state'].isin(['tpbase', 'fn'])].copy()
base['base_id'] = base['MatchId'].apply(lambda x: x[0])
base['comp_id'] = base['MatchId'].apply(lambda x: x[1])

# Separate out the variants from the comparison VCF and add new columns of the base/comp ids
comp = data[data['state'].isin(['tp', 'fp'])].copy()
comp['base_id'] = comp['MatchId'].apply(lambda x: x[0])
comp['comp_id'] = comp['MatchId'].apply(lambda x: x[1])

# Merge the base/comparison variants
combined = pd.merge(base, comp, left_on='base_id', right_on='comp_id', suffixes=('_base', '_comp'))

# How many comp variants matched to multiple base variants?
counts1 = combined['base_id_comp'].value_counts()
print('multi-matched comp count', (counts1 != 1).sum())

# How many base variants matched to multiple comp variants?
counts2 = combined['comp_id_base'].value_counts()
print('multi-matched base count', (counts2 != 1).sum())

# For every tp-base, I can check its GTMatch to figure out concordance/discordance
#.. but that doesn't solve the ALstate lookups..
# So I need to use a vcf2df 
# and for every fn, I add that to discorant
# and for every fp, I add that to discordant
# and then I can add to concordant every reference homozygous base_vcf
# And then I can add to missing every Missing
