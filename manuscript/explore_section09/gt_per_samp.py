import sys
import truvari
import pysam
import numpy as np

fh = pysam.VariantFile(sys.argv[1])
nsamp = len(fh.header.samples)
cnts = np.zeros((len(truvari.GT), nsamp))

for entry in fh:
    if entry.chrom in ['chrX', 'chrY']:
        continue
    for i, samp in enumerate(entry.samples.values()):
        cnts[truvari.get_gt(samp['GT']).value, i] += 1

import joblib
joblib.dump(cnts, 'gtcnts.jl')
