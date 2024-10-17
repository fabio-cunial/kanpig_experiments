import joblib
import pandas as pd
import json
from scipy.stats import pearsonr
import sys
import truvari

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

def process(df):
    cor = get_pairs(df)
    x = df['state'].value_counts().to_dict()
    perf = truvari.performance_metrics(**x)
    result = dict(zip(['precision', 'recall', 'f1'], perf))
    result['correlation'] = cor.statistic
    result['pvalue'] = cor.pvalue
    return result


results = {}
results["ALL"] = process(orig)

view = orig[orig['svtype'] == 'DEL'].copy()
results["DEL"] = process(view)

view = orig[orig['svtype'] == 'INS'].copy()
results["INS"] = process(view)
print(json.dumps(results, indent=4))


