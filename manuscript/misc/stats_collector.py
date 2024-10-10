import sys
import joblib
import pandas as pd
import pysam

parts = []
meta = pd.read_csv("../metadata.txt", sep='\t')
meta.set_index('Sample', inplace=True)
print("Pipeline\tSample\tPopulation\tSuper-Population\tVariant Count\tHet/Hom Ratio")
for program in ['1kgp', 'kanpig']:
    for population in ['AFR', 'AMR', 'EAS', 'EUR', 'SAS']:
        vcf_fn = f"{program}.{population}.vcf.gz"
        stats = joblib.load(f"stats.{vcf_fn}.jl")
        vcf = pysam.VariantFile(vcf_fn)
        for idx, sample in enumerate(vcf.header.samples):
            pop = meta.loc[sample]['Population']
            print(program, sample, population, pop, stats[idx].sum(), stats[idx][0] / stats[idx][1], sep='\t')
