import joblib
import pandas as pd
from scipy.stats import pearsonr
import sys

def get_af(x):
    if isinstance(x, tuple):
        x = x[0]
    return x

def get_pairs(data):
    #data.reset_index(inplace=True)
    # Separate out the variants from the base VCF and add new columns of the base/comp ids
    data['AF'] = data['AF'].apply(get_af)

    base = data[data['state'].isin(['tpbase'])].copy()
    base['base_id'] = base['MatchId'].apply(lambda x: x[0])
    base['comp_id'] = base['MatchId'].apply(lambda x: x[1])

    # Separate out the variants from the comparison VCF and add new columns of the base/comp ids
    comp = data[data['state'].isin(['tp'])].copy()
    comp['base_id'] = comp['MatchId'].apply(lambda x: x[0])
    comp['comp_id'] = comp['MatchId'].apply(lambda x: x[1])

    # Merge the base/comparison variants
    combined = pd.merge(base, comp, left_on='base_id', right_on='comp_id', suffixes=('_base', '_comp'))
    return pearsonr(combined['AF_base'], combined['AF_comp'])

orig = joblib.load(sys.argv[1])
print(get_pairs(orig))
#for key in orig:
    #print(key, get_pairs(orig[key]))


