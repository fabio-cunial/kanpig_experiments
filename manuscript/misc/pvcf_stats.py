import sys
import numpy as np
import pysam
import truvari
import joblib
out_name = sys.argv[2]
vcf = pysam.VariantFile(sys.argv[1])

counts = np.zeros((len(vcf.header.samples), 2))

for entry in vcf:
    if entry.chrom in ['chrX', 'chrY', 'chrM']:
        continue
    for idx, sample in enumerate(entry.samples.values()):
        gt = truvari.get_gt(sample['GT'])
        if gt == truvari.GT.HOM:
            g = 1
        elif gt == truvari.GT.HET:
            g = 0
        else:
            continue
        counts[idx][g] += 1

joblib.dump(counts, out_name)
